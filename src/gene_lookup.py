"""Reference-based gene -> transcript -> exon -> coding codon workflow.

Gene symbols and Ensembl IDs are resolved against Ensembl's CURRENT REST
reference database, not a frozen genome version. GRCh38 is required for this
human-only workflow. Translation and exon mapping are checked against the
returned CDS length before any genotype-assay design may be attempted.

A genome-wide primer-specificity search and clinical genotype calls are NOT
performed. No user-supplied patient sequences are sent to Ensembl.
"""
from collections import defaultdict
from urllib.parse import quote, urlencode
from urllib.request import build_opener, HTTPSHandler, HTTPRedirectHandler, Request
from urllib.error import HTTPError, URLError
import json
import re

from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.Data.CodonTable import TranslationError
from .sequence import load_sequence

ENSEMBL = "https://rest.ensembl.org"
SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,39}$")
GENE_ID_RE = re.compile(r"^ENSG[0-9]{11}$", re.I)
TX_ID_RE = re.compile(r"^ENST[0-9]{11}$", re.I)
CHROM_RE = re.compile(r"^(?:[0-9]{1,2}|X|Y)$", re.I)
MAX_RESPONSE = 8 * 1024 * 1024
MAX_CDS = 150_000
MAX_GENE_TRANSCRIPTS = 1500


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError("Unexpected Ensembl redirect blocked.")


_opener = build_opener(HTTPSHandler(), _NoRedirect())


def _get_json(resource, params=None):
    """Fixed-host HTTPS GET; no user-supplied URL, no redirects, bounded response."""
    if not resource.startswith("/") or "//" in resource or ".." in resource:
        raise ValueError("Invalid Ensembl resource path.")
    uri = ENSEMBL + resource
    if params:
        uri += "?" + urlencode(params)
    request = Request(uri, headers={"Accept": "application/json", "User-Agent": "Genotype-PCR-RFLP-research/0.3"})
    try:
        with _opener.open(request, timeout=12) as response:
            data = response.read(MAX_RESPONSE + 1)
            if len(data) > MAX_RESPONSE:
                raise ValueError("Ensembl response exceeds supported size.")
            return json.loads(data.decode("utf-8"))
    except HTTPError as exc:
        if exc.code == 404:
            raise ValueError("Gene/transcript/sequence was not found in Ensembl.") from exc
        raise ValueError(f"Ensembl returned HTTP {exc.code}. Please retry later.") from exc
    except (TimeoutError, URLError, OSError) as exc:
        raise ValueError("Ensembl is unavailable or the request timed out. Internet access is required.") from exc
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("Ensembl returned an unreadable reference response.") from exc


def search_gene(identifier):
    """Lookup one human gene by symbol or Ensembl gene ID."""
    name = str(identifier).strip()
    if GENE_ID_RE.fullmatch(name):
        endpoint = "/lookup/id/" + quote(name.upper(), safe="")
    elif SYMBOL_RE.fullmatch(name):
        endpoint = "/lookup/symbol/homo_sapiens/" + quote(name, safe="")
    else:
        raise ValueError("Enter a human gene symbol or Ensembl gene ID (ENSG...).")
    data = _get_json(endpoint, {"expand": 1})
    if data.get("object_type") != "Gene" or data.get("species") != "homo_sapiens":
        raise ValueError("The requested identifier did not resolve to a human gene.")
    if data.get("assembly_name") not in (None, "GRCh38"):
        raise ValueError("This workflow requires GRCh38 reference annotation.")
    txs = data.get("Transcript", [])
    if not isinstance(txs, list) or len(txs) > MAX_GENE_TRANSCRIPTS:
        raise ValueError("The gene has too many transcripts for this interactive lookup.")
    transcripts = []
    for t in txs:
        if not t.get("Translation") or t.get("biotype") != "protein_coding":
            continue
        if not TX_ID_RE.fullmatch(str(t.get("id", ""))):
            continue
        if not t.get("Exon") or str(t.get("seq_region_name", "")).upper() == "MT":
            continue
        transcripts.append({
            "id": t["id"], "version": t.get("version"), "label": t.get("display_name", t["id"]),
            "is_canonical": bool(t.get("is_canonical", 0)),
            "coding_length_aa": t.get("Translation", {}).get("length"),
            "exon_count": len(t["Exon"]),
            "strand": t["strand"],
            "seq_region_name": t["seq_region_name"],
            "assembly_name": t.get("assembly_name", data.get("assembly_name")),
            "annotation": t,
        })
    transcripts.sort(key=lambda t: (not t["is_canonical"], -(t["coding_length_aa"] or 0), t["id"]))
    if not transcripts:
        raise ValueError("No protein-coding transcript with exon and translation annotations was found.")
    return {
        "gene_id": data["id"],
        "gene_symbol": data.get("display_name", name),
        "assembly": "GRCh38 (Ensembl default reference)",
        "chromosome": data.get("seq_region_name"),
        "transcripts": transcripts,
        "reference_provider": "Ensembl REST",
        "reference_url": ENSEMBL + endpoint,
    }


def get_cds(transcript_id):
    if not TX_ID_RE.fullmatch(str(transcript_id)):
        raise ValueError("Select a valid Ensembl ENST transcript.")
    response = _get_json("/sequence/id/" + transcript_id, {"type": "cds"})
    sequence = load_sequence(str(response.get("seq", "")))
    if len(sequence) > MAX_CDS or len(sequence) % 3:
        raise ValueError("The CDS length is unsupported or not divisible by three.")
    return sequence


def coding_map(transcript, cds):
    """Map each 1-based CDS nucleotide to its genomic coordinate and exon number.

    Translation.start/end are Ensembl genomic coordinates of CDS boundaries,
    and Exon coordinates are included in the expanded transcript lookup.
    On negative-strand transcripts genomic coordinates descend in CDS order.
    Fails closed if exon spans cannot reconcile to the supplied CDS.
    """
    strand = int(transcript.get("strand", 0))
    if strand not in (1, -1):
        raise ValueError("Transcript strand is invalid.")
    translation = transcript.get("Translation") or {}
    try:
        first = int(translation["start"])
        last = int(translation["end"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Selected transcript has no genomic translation bounds.") from exc
    low, high = sorted((first, last))
    exons = sorted(transcript.get("Exon", []), key=lambda exon: int(exon["start"]), reverse=strand == -1)
    if not exons:
        raise ValueError("No exons available for this transcript.")
    coords = []
    exon_labels = []
    for exon_index, exon in enumerate(exons, 1):
        e_start, e_end = int(exon["start"]), int(exon["end"])
        s, e = max(e_start, low), min(e_end, high)
        if s > e:
            continue
        stretch = range(s, e + 1) if strand == 1 else range(e, s - 1, -1)
        for genomic in stretch:
            coords.append(genomic)
            exon_labels.append(exon_index)
        if len(coords) > MAX_CDS:
            raise ValueError("Coding mapping exceeds safe length.")
    if len(coords) != len(cds):
        raise ValueError(
            f"Ensembl exon mapping spans {len(coords)} bases but returned CDS has {len(cds)}. "
            "Cannot assign safe codon coordinates for this transcript."
        )
    return coords, exon_labels


def codon_table(transcript, cds, max_rows=None):
    coords, exon_labels = coding_map(transcript, cds)
    rows = []
    for offset in range(0, len(cds), 3):
        codon = cds[offset:offset + 3]
        aa = str(Seq(codon).translate(table=1))
        rows.append({
            "codon_number": offset // 3 + 1,
            "reference_codon": codon,
            "reference_amino_acid": aa,
            "exon_numbers": sorted(set(exon_labels[offset:offset + 3])),
            "genomic_positions": coords[offset:offset + 3],
            "cds_positions_1based": [offset + 1, offset + 2, offset + 3],
            "splice_junction_codon": len(set(exon_labels[offset:offset + 3])) > 1,
        })
    if max_rows is not None:
        return rows[:max_rows]
    return rows


def alternative_codons(codon):
    """All nine single-nucleotide alternate codons with predicted amino acid."""
    codon = str(codon).upper()
    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        raise ValueError("Expected a three-base reference coding codon.")
    options = []
    reference_aa = str(Seq(codon).translate(table=1))
    for i in range(3):
        for alt in "ACGT":
            if alt == codon[i]:
                continue
            other = codon[:i] + alt + codon[i + 1:]
            aa = str(Seq(other).translate(table=1))
            options.append({
                "mutant_codon": other, "codon_base": i + 1,
                "reference_base": codon[i], "alternate_base": alt,
                "mutant_amino_acid": aa,
                "protein_effect": "synonymous" if aa == reference_aa else "stop gained" if aa == "*" else "missense/nonsense change",
            })
    return options


def fetch_genomic_assay(transcript, cds, codon_number, mutant_codon, flank_bp=300):
    """Return a genomic DNA reference flank plus local SNV for current engine.

    Ensembl supplies genomic plus strand. Negative-strand transcripts are
    reverse-complemented so the local mutation has the same allele as CDS.
    """
    if int(transcript.get("strand", 0)) not in (1, -1):
        raise ValueError("Invalid selected transcript strand.")
    if not isinstance(codon_number, int) or codon_number < 1 or codon_number * 3 > len(cds):
        raise ValueError("Codon position is outside the selected transcript's CDS.")
    if flank_bp < 200 or flank_bp > 1000:
        raise ValueError("Flank must be 200–1000 bp.")
    codon = cds[(codon_number - 1)*3:codon_number*3]
    options = alternative_codons(codon)
    change = next((v for v in options if v["mutant_codon"] == mutant_codon.upper()), None)
    if not change:
        raise ValueError("Alternate codon must differ from reference by exactly one base.")
    map_coords, exon_labels = coding_map(transcript, cds)
    idx = 3 * (codon_number - 1) + change["codon_base"] - 1
    genomic_pos = map_coords[idx]
    chrom = str(transcript.get("seq_region_name", ""))
    if not CHROM_RE.fullmatch(chrom):
        raise ValueError("Unsupported chromosome identifier for human genomic region lookup.")
    lo, hi = max(1, genomic_pos - flank_bp), genomic_pos + flank_bp
    region = f"{chrom}:{lo}..{hi}:1"
    resp = _get_json("/sequence/region/human/" + quote(region, safe=":."), {})
    plus_sequence = load_sequence(str(resp.get("seq", "")))
    if len(plus_sequence) != hi - lo + 1:
        raise ValueError("Ensembl returned a region of an unexpected length.")
    negative = transcript["strand"] == -1
    dna = str(Seq(plus_sequence).reverse_complement()) if negative else plus_sequence
    local_index = hi - genomic_pos if negative else genomic_pos - lo
    ref, alt = change["reference_base"], change["alternate_base"]
    if dna[local_index] != ref:
        raise ValueError("CDS/reference genomic flank disagree at the target base; refusing assay design.")
    alt_seq = dna[:local_index] + alt + dna[local_index + 1:]
    aa_before = str(Seq(codon).translate(table=1))
    aa_after = change["mutant_amino_acid"]
    metadata = {
        "provider": "Ensembl REST; current reference release",
        "assembly": "GRCh38",
        "gene": transcript.get("gene_symbol", ""),
        "transcript_id": transcript["id"],
        "transcript_version": transcript.get("version"),
        "exon_number": exon_labels[idx],
        "codon_number": codon_number,
        "coding_position_c": idx + 1,
        "reference_codon": codon,
        "alternate_codon": mutant_codon.upper(),
        "reference_aa": aa_before,
        "alternate_aa": aa_after,
        "genomic_chromosome": chrom,
        "genomic_position": genomic_pos,
        "genomic_interval": f"{chrom}:{lo}-{hi}",
        "transcript_strand": transcript["strand"],
        "local_mutation": f"{ref}{local_index+1}{alt}",
        "reference_sequence_source": f"{ENSEMBL}/sequence/region/human/{region}",
        "mapping_status": "CDS-to-EXON-CHECKED; GENOMIC-BASE-CHECKED",
        "limitation": "A predicted reference-to-mutant change, not a genotype call from an actual sample. Primer off-target specificity and experimental verification are not performed.",
    }
    return {"sequence": dna, "mutation": metadata["local_mutation"], "metadata": metadata}
