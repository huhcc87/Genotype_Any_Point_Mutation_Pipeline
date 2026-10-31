"""Independent UCSC RefSeq gene -> exon -> coding DNA reference workflow.

UCSC ncbiRefSeqCurated genePred track provides accession-versioned NM_
transcripts; coordinates are UCSC zero-based half-open genomic intervals.
This module uses assembly-matched hg38/hg19 sequence and converts boundaries
to 1-based inclusive for the shared codon mapper. It never merges Ensembl
ENST and NCBI RefSeq NM_ transcript coordinate systems.
"""
import json
import re
import time
from urllib.parse import urlencode
from urllib.request import Request, build_opener, HTTPSHandler, HTTPRedirectHandler
from urllib.error import HTTPError, URLError

from Bio.Seq import Seq
from .sequence import load_sequence

_UCSC = "https://api.genome.ucsc.edu"
_MYGENE = "https://mygene.info"
_GENOMES = {"GRCh38": "hg38", "GRCh37": "hg19"}
_SYMBOL = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,39}$")
_ENSG = re.compile(r"^ENSG[0-9]{11}$", re.I)
_NM = re.compile(r"^NM_[0-9]+\.[0-9]+$")
_CHROM = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
_MAX_BYTES = 5_000_000
_MAX_SPAN = 700_000
_MAX_TRANSCRIPTS = 300


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError("Unexpected reference-service redirect.")


_opener = build_opener(HTTPSHandler(), _NoRedirect())


def _read_json(host, path, params):
    if host not in (_UCSC, _MYGENE) or not path.startswith("/") or ".." in path:
        raise ValueError("Invalid fixed reference-service endpoint.")
    url = host + path + "?" + urlencode(params)
    req = Request(url, headers={"Accept": "application/json", "User-Agent": "Genotype-Research-Pipeline/0.5"})
    for attempt in range(2):
        try:
            with _opener.open(req, timeout=15) as response:
                body = response.read(_MAX_BYTES + 1)
                if len(body) > _MAX_BYTES:
                    raise ValueError("Reference service response exceeded the safe limit.")
                value = json.loads(body.decode("utf-8"))
                if not isinstance(value, dict):
                    raise ValueError("Expected a JSON object from the reference service.")
                if value.get("error"):
                    raise ValueError(f"Reference service error for {path}: {str(value['error'])[:300]}")
                return value
        except HTTPError as exc:
            if exc.code in (429, 500, 502, 503, 504) and attempt == 0:
                time.sleep(0.5)
                continue
            raise ValueError(f"{host}{path} returned HTTP {exc.code}") from exc
        except (OSError, TimeoutError, URLError) as exc:
            if attempt == 0:
                time.sleep(0.5)
                continue
            raise ValueError(f"{host}{path} timed out or was unreachable: {exc}") from exc
        except (UnicodeError, json.JSONDecodeError) as exc:
            raise ValueError(f"{host}{path} returned invalid JSON.") from exc
    raise ValueError(f"Reference service failed for {host}{path}.")


def _resolve_gene_locus(symbol, assembly):
    """Look up exact human gene symbol and selected assembly's locus via MyGene."""
    q = "ensembl.gene:" + symbol if _ENSG.fullmatch(symbol) else "symbol:" + symbol
    data = _read_json(_MYGENE, "/v3/query", {
        "q": q, "species": "9606", "size": 10,
        "fields": "symbol,taxid,ensembl,genomic_pos,genomic_pos_hg19"
    })
    for hit in data.get("hits", []):
        if not isinstance(hit, dict) or hit.get("taxid") != 9606:
            continue
        if not _ENSG.fullmatch(symbol) and str(hit.get("symbol", "")).upper() != symbol.upper():
            continue
        if _ENSG.fullmatch(symbol):
            ens = hit.get("ensembl", [])
            ens = [ens] if isinstance(ens, dict) else ens if isinstance(ens, list) else []
            if not any(str(x.get("gene","")).split(".")[0].upper() == symbol.upper() for x in ens):
                continue
        genome_field = "genomic_pos" if assembly == "GRCh38" else "genomic_pos_hg19"
        positions = hit.get(genome_field, [])
        positions = [positions] if isinstance(positions, dict) else positions if isinstance(positions, list) else []
        for locus in positions:
            if not isinstance(locus, dict):
                continue
            chrom = "chr" + str(locus.get("chr","")).removeprefix("chr")
            try:
                start, end = int(locus["start"]), int(locus["end"])
            except (KeyError, ValueError, TypeError):
                continue
            if _CHROM.fullmatch(chrom) and 1 <= start <= end and end-start+1 <= _MAX_SPAN:
                return {"gene_symbol": hit.get("symbol",symbol),
                        "chromosome": chrom, "start_1based": start,
                        "end_1based": end, "locus_source": "MyGene.info"}
    raise ValueError("No exact human gene locus found for this assembly in MyGene.info.")


def _parse_bed_positions(value, count, label):
    if isinstance(value, str):
        positions = [int(n) for n in value.strip(",").split(",") if n]
    elif isinstance(value, list):
        positions = [int(n) for n in value]
    else:
        raise ValueError(f"Missing {label} exon boundaries.")
    if len(positions) != count:
        raise ValueError(f"RefSeq {label} exon coordinates do not match exon count.")
    return positions


def _normalize_refseq(row, expected_symbol, expected_chrom):
    """Strict RefSeq genePred coordinate validation and normalized transcript."""
    if not isinstance(row, dict):
        return None
    if str(row.get("name2","")).upper() != expected_symbol.upper():
        return None
    name = str(row.get("name",""))
    if not _NM.fullmatch(name) or row.get("chrom") != expected_chrom:
        return None
    if row.get("strand") not in ("+", "-"):
        return None
    if row.get("cdsStartStat") not in (None, "cmpl") or row.get("cdsEndStat") not in (None, "cmpl"):
        return None
    try:
        exons = int(row["exonCount"])
        start, end = int(row["txStart"]), int(row["txEnd"])
        cds_start, cds_end = int(row["cdsStart"]), int(row["cdsEnd"])
        starts = _parse_bed_positions(row["exonStarts"], exons, "start")
        ends = _parse_bed_positions(row["exonEnds"], exons, "end")
    except (KeyError, ValueError, TypeError):
        return None
    if not (0 <= start < cds_start < cds_end <= end and exons > 0 and exons < 500):
        # Some valid first-codon transcripts have cdsStart==txStart.
        if not (0 <= start <= cds_start < cds_end <= end and 0 < exons < 500):
            return None
    if any(a < start or b > end or a >= b for a,b in zip(starts,ends)):
        return None
    if any(starts[i] >= starts[i+1] or ends[i] > starts[i+1] for i in range(exons-1)):
        return None
    strand = 1 if row["strand"] == "+" else -1
    ensembl_like = {
        "id": name, "strand": strand,
        "seq_region_name": expected_chrom[3:],
        "biotype": "protein_coding",
        "source": "UCSC RefSeq Curated",
        "Exon": [{"start": a+1, "end": z} for a,z in zip(starts,ends)],
        "Translation": {"start": cds_start+1, "end": cds_end},
        "refseq_track": {"chrom":expected_chrom,"txStart":start,"txEnd":end,
                         "cdsStart":cds_start,"cdsEnd":cds_end},
    }
    return {
        "id":name, "version":int(name.split(".")[-1]), "label":name,
        "is_canonical":False, "coding_length_aa":None, "exon_count":exons,
        "strand":strand, "seq_region_name":expected_chrom[3:],
        "assembly_name":"Reference annotation must be checked",
        "annotation": ensembl_like,
        "provider":"UCSC RefSeq Curated",
    }


def search_refseq_gene(gene_symbol, assembly="GRCh38"):
    """Resolve gene locus then query curated NM_ transcripts on same assembly."""
    name = str(gene_symbol).strip()
    if assembly not in _GENOMES or not (_SYMBOL.fullmatch(name) or _ENSG.fullmatch(name)):
        raise ValueError("Enter a human gene symbol/ENSG and choose GRCh38 or GRCh37.")
    locus = _resolve_gene_locus(name, assembly)
    start0 = max(0, locus["start_1based"] - 501)
    end0 = locus["end_1based"] + 500
    if end0 - start0 > _MAX_SPAN + 1000:
        raise ValueError("Gene locus is larger than supported UCSC reference interval.")
    result = _read_json(_UCSC, "/getData/track", {
        "genome": _GENOMES[assembly], "track": "ncbiRefSeqCurated",
        "chrom":locus["chromosome"],"start":start0,"end":end0,
        "maxItemsOutput": 1500,
    })
    rows = result.get("ncbiRefSeqCurated")
    if not isinstance(rows, list):
        raise ValueError("UCSC returned no curated RefSeq track records.")
    transcripts = []
    for row in rows:
        tx = _normalize_refseq(row, locus["gene_symbol"], locus["chromosome"])
        if tx is not None:
            transcripts.append(tx)
        if len(transcripts) >= _MAX_TRANSCRIPTS:
            break
    if not transcripts:
        raise ValueError(
            "No matching curated NM_ transcript for this gene/assembly. "
            "Try Ensembl mode or use a validated reference FASTA."
        )
    transcripts.sort(key=lambda tx: (-tx["exon_count"], tx["id"]))
    return {
        "gene_id": name if _ENSG.fullmatch(name) else locus["gene_symbol"],
        "gene_locus":{"start_1based":locus["start_1based"],"end_1based":locus["end_1based"]},
        "gene_locus_source":"MyGene.info genomic position (selected assembly)",
        "gene_symbol": locus["gene_symbol"],
        "assembly":assembly,"chromosome":locus["chromosome"],
        "transcripts":transcripts,
        "reference_provider":"UCSC Genome Browser ncbiRefSeqCurated / MyGene locus",
        "reference_url":_UCSC + "/getData/track",
        "query":name,
        "lookup_status":"UCSC_REFSEQ_CURATED",
        "lookup_warning":(
            "Transcript models are NCBI RefSeq NM_ accessions, not Ensembl ENST. "
            "Exon and codon numbering is specific to selected RefSeq NM_ transcript. "
            "No identity between RefSeq and Ensembl transcripts is assumed."
        )
    }


def get_refseq_cds(transcript, assembly="GRCh38"):
    """Fetch genome DNA for the selected NM_ and reconstruct its spliced CDS."""
    if assembly not in _GENOMES or transcript.get("source") != "UCSC RefSeq Curated":
        raise ValueError("Expected UCSC RefSeq transcript annotation with selected assembly.")
    track = transcript["refseq_track"]
    sequence = fetch_refseq_region(track["chrom"],track["txStart"],track["txEnd"],assembly)
    exon_chunks = []
    for exon in transcript["Exon"]:
        left = max(exon["start"] - 1, track["cdsStart"])
        right = min(exon["end"], track["cdsEnd"])
        if left < right:
            exon_chunks.append(sequence[left-track["txStart"]:right-track["txStart"]])
    cds_plus = "".join(exon_chunks)
    cds = str(Seq(cds_plus).reverse_complement()) if transcript["strand"] == -1 else cds_plus
    if not cds or len(cds) % 3:
        raise ValueError("Reference CDS is not divisible by 3; transcript cannot be used.")
    if len(cds) > 150_000:
        raise ValueError("RefSeq CDS exceeds supported sequence length.")
    # Do not silently accept a RefSeq model with clearly inconsistent translation.
    if not cds.startswith("ATG"):
        raise ValueError("RefSeq CDS lacks an initial ATG; transcript may be incomplete or use alternate translation.")
    return cds


def fetch_refseq_region(chrom, start0, end0, assembly="GRCh38"):
    if assembly not in _GENOMES or not _CHROM.fullmatch(chrom):
        raise ValueError("Unsupported assembly or chromosome.")
    if isinstance(start0,bool) or isinstance(end0,bool) or not (0 <= start0 < end0 <= 1_000_000_000):
        raise ValueError("Invalid zero-based genomic interval.")
    if end0-start0 > _MAX_SPAN:
        raise ValueError("Reference sequence interval exceeds safe size.")
    response = _read_json(_UCSC, "/getData/sequence", {
        "genome":_GENOMES[assembly],"chrom":chrom,"start":start0,"end":end0,
    })
    dna = load_sequence(str(response.get("dna","")))
    if len(dna) != end0-start0:
        raise ValueError("UCSC returned incorrect genome sequence interval length.")
    return dna


def fetch_refseq_assay(transcript, cds, codon_number, alternate_codon,
                       flank_bp=300, assembly="GRCh38"):
    """Use a validated codon map and UCSC genomic reference to create local SNV."""
    from .gene_lookup import coding_map, alternative_codons
    if transcript.get("source") != "UCSC RefSeq Curated":
        raise ValueError("Expected UCSC RefSeq transcript.")
    if not isinstance(codon_number,int) or not (1 <= codon_number <= len(cds)//3):
        raise ValueError("Codon is outside coding sequence.")
    codon = cds[(codon_number-1)*3:codon_number*3]
    variant = next((v for v in alternative_codons(codon)
                   if v["mutant_codon"] == alternate_codon),None)
    if variant is None:
        raise ValueError("Mutation must change exactly one base of reference codon.")
    coordinates, exon_labels = coding_map(transcript,cds)
    idx = (codon_number-1)*3 + variant["codon_base"]-1
    chrom = "chr"+str(transcript["seq_region_name"])
    position = coordinates[idx]
    left0=max(0, position-flank_bp-1)
    right0=position+flank_bp
    plus=fetch_refseq_region(chrom,left0,right0,assembly)
    reverse=transcript["strand"] == -1
    sequence=str(Seq(plus).reverse_complement()) if reverse else plus
    local=right0-position if reverse else position-1-left0
    ref,alt=variant["reference_base"],variant["alternate_base"]
    if sequence[local]!=ref:
        raise ValueError("RefSeq CDS and UCSC genomic DNA disagree at variant: assay blocked.")
    return {
        "sequence":sequence,
        "mutation":f"{ref}{local+1}{alt}",
        "metadata":{
            "provider":"UCSC Genome Browser RefSeq Curated",
            "assembly":assembly,"transcript_id":transcript["id"],
            "transcript_version":transcript["id"].split(".")[-1],
            "exon_number":exon_labels[idx],"codon_number":codon_number,
            "coding_position_c":idx+1,"reference_codon":codon,
            "alternate_codon":alternate_codon,
            "reference_aa":str(Seq(codon).translate(table=1)),
            "alternate_aa":variant["mutant_amino_acid"],
            "genomic_chromosome":chrom,"genomic_position":position,
            "genomic_interval":f"{chrom}:{left0+1}-{right0}",
            "transcript_strand":transcript["strand"],
            "local_mutation":f"{ref}{local+1}{alt}",
            "mapping_status":"REFSEQ_CDS_AND_GENOMIC_BASE_VERIFIED",
            "limitation":"Research-use reference prediction only; genome-wide primer specificity and experimental genotype calling remain unverified.",
        },
    }
