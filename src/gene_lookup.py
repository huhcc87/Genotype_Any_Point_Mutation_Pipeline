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
HOSTS = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
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


def _get_json(resource, params=None, assembly="GRCh38"):
    """Bounded Ensembl HTTPS request with short retries for transient failures."""
    import time
    if not resource.startswith("/") or "//" in resource or ".." in resource:
        raise ValueError("Invalid Ensembl resource path.")
    if assembly not in HOSTS:
        raise ValueError("Supported reference assemblies are GRCh38 and GRCh37.")
    uri = HOSTS[assembly] + resource
    if params:
        uri += "?" + urlencode(params)
    request = Request(uri, headers={
        "Accept": "application/json",
        "Content-Type": "application/json",
        "User-Agent": "Genotype-PCR-RFLP-research/0.4",
    })
    last_error = None
    for attempt in range(2):
        try:
            with _opener.open(request, timeout=16) as response:
                data = response.read(MAX_RESPONSE + 1)
                if len(data) > MAX_RESPONSE:
                    raise ValueError("Ensembl response too large; choose a more specific transcript.")
                return json.loads(data.decode("utf-8"))
        except HTTPError as exc:
            if exc.code == 404:
                raise ValueError("Gene, transcript or reference sequence was not found in Ensembl.") from exc
            if exc.code in (429, 500, 502, 503, 504) and attempt == 0:
                last_error = f"HTTP {exc.code}"
                time.sleep(1)
                continue
            raise ValueError(
                f"Ensembl returned HTTP {exc.code}. Check the Ensembl service status, or retry later."
            ) from exc
        except (TimeoutError, URLError, OSError) as exc:
            last_error = str(getattr(exc, "reason", exc))[:160]
            if attempt == 0:
                time.sleep(1)
                continue
            raise ValueError(
                "Cannot reach Ensembl REST. Check internet, VPN/proxy/firewall or service availability. "
                f"Network detail: {last_error}. You may still use Paste sequence or WT/Mutant input offline."
            ) from exc
        except (UnicodeError, json.JSONDecodeError) as exc:
            raise ValueError("Ensembl returned an unreadable reference response.") from exc
    raise ValueError(f"Reference service unavailable: {last_error or 'unknown reason'}")


def search_gene(identifier, assembly="GRCh38"):
    """Fast lookup: one gene record plus lightweight coding transcript features.

    Do not use expand=1 on a gene: very large genes such as TTN can time out
    or return multi-megabyte nested exon data.
    """
    name = str(identifier).strip()
    if GENE_ID_RE.fullmatch(name):
        endpoint = "/lookup/id/" + quote(name.upper(), safe="")
    elif SYMBOL_RE.fullmatch(name):
        endpoint = "/lookup/symbol/homo_sapiens/" + quote(name, safe="")
    else:
        raise ValueError("Enter a human gene symbol or Ensembl ENSG gene ID.")
    gene = _get_json(endpoint, assembly=assembly)
    if not isinstance(gene, dict) or gene.get("object_type") != "Gene" or gene.get("species") != "homo_sapiens":
        raise ValueError("Identifier did not resolve to a human gene.")
    if gene.get("assembly_name") not in (None, assembly):
        raise ValueError("Expected the GRCh38 human reference.")
    gene_id = str(gene.get("id", ""))
    if not GENE_ID_RE.fullmatch(gene_id):
        raise ValueError("Ensembl did not return a valid human gene ID.")

    features = _get_json(
        "/overlap/id/" + quote(gene_id, safe=""),
        {"feature": "transcript", "biotype": "protein_coding"}, assembly=assembly,
    )
    if not isinstance(features, list) or len(features) > MAX_GENE_TRANSCRIPTS:
        raise ValueError("Ensembl transcript response is invalid or exceeds the supported limit.")
    canonical = str(gene.get("canonical_transcript") or "").split(".")[0]
    transcripts = []
    for t in features:
        if not isinstance(t, dict):
            continue
        # overlap/id returns overlaps from neighbouring genes too.
        parent = str(t.get("Parent") or t.get("gene_id") or "")
        if parent != gene_id:
            continue
        tid = str(t.get("id", ""))
        if not TX_ID_RE.fullmatch(tid):
            continue
        if t.get("biotype") != "protein_coding":
            continue
        if str(t.get("seq_region_name")) == "MT":
            continue
        transcripts.append({
            "id": tid,
            "version": t.get("version"),
            "label": t.get("external_name") or t.get("display_name") or tid,
            "is_canonical": tid == canonical,
            "coding_length_aa": None,
            "exon_count": None,
            "strand": t.get("strand"),
            "seq_region_name": t.get("seq_region_name"),
            "assembly_name": t.get("assembly_name", assembly),
            # Detail is intentionally loaded only after transcript selection.
            "annotation": None,
        })
    transcripts.sort(key=lambda t: (not t["is_canonical"], t["id"]))
    if not transcripts:
        raise ValueError(
            "No matching protein-coding transcripts were returned. "
            "Check the gene identifier or try another reference sequence."
        )
    return {
        "gene_id": gene_id,
        "gene_symbol": gene.get("display_name", name),
        "assembly": assembly,
        "chromosome": gene.get("seq_region_name"),
        "transcripts": transcripts,
        "reference_provider": "Ensembl REST",
        "reference_url": HOSTS[assembly] + endpoint,
        "query": name,
    }


def get_transcript_details(transcript_id, expected_gene_id=None, assembly="GRCh38"):
    """Expand one ENST only; validate gene membership and coding annotations."""
    tid = str(transcript_id).strip().upper()
    if not TX_ID_RE.fullmatch(tid):
        raise ValueError("Select a valid ENST transcript ID.")
    transcript = _get_json("/lookup/id/" + quote(tid, safe=""), {"expand": 1}, assembly=assembly)
    if not isinstance(transcript, dict) or transcript.get("object_type") != "Transcript":
        raise ValueError("Ensembl returned no transcript detail.")
    if transcript.get("id") != tid or transcript.get("biotype") != "protein_coding":
        raise ValueError("The selected transcript is not protein coding.")
    if expected_gene_id and transcript.get("Parent") != expected_gene_id:
        raise ValueError("Selected transcript does not belong to the requested gene.")
    if transcript.get("assembly_name") not in (None, assembly):
        raise ValueError("Transcript does not use GRCh38.")
    if not transcript.get("Translation") or not transcript.get("Exon"):
        raise ValueError("This transcript has no usable CDS translation or exon annotations.")
    if transcript.get("seq_region_name") == "MT":
        raise ValueError("Mitochondrial transcripts require a separate codon table and are unsupported.")
    return transcript


def get_cds(transcript_id, assembly="GRCh38"):
    if not TX_ID_RE.fullmatch(str(transcript_id)):
        raise ValueError("Select a valid Ensembl ENST transcript.")
    response = _get_json("/sequence/id/" + transcript_id, {"type": "cds"}, assembly=assembly)
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


def fetch_genomic_assay(transcript, cds, codon_number, mutant_codon, flank_bp=300, assembly="GRCh38"):
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
    resp = _get_json("/sequence/region/human/" + quote(region, safe=":."), {}, assembly=assembly)
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
        "assembly": assembly,
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
        "reference_sequence_source": f"{HOSTS[assembly]}/sequence/region/human/{region}",
        "mapping_status": "CDS-to-EXON-CHECKED; GENOMIC-BASE-CHECKED",
        "limitation": "A predicted reference-to-mutant change, not a genotype call from an actual sample. Primer off-target specificity and experimental verification are not performed.",
    }
    return {"sequence": dna, "mutation": metadata["local_mutation"], "metadata": metadata}
