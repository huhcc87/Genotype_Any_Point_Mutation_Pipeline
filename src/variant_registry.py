"""Research-only NCBI variant record lookup for a *genomic* SNV.

A returned rsID is an allele-associated identifier, not a clinical genotype.
ClinVar esearch results are candidate records until checked for exact allele
identity; an empty database response is not evidence that a variant is novel.
Only public coordinates/alleles are sent to NCBI, never sample sequences.
"""
import json
import re
from datetime import datetime, timezone
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, build_opener, HTTPSHandler, HTTPRedirectHandler
from Bio.Seq import Seq

NCBI_VARIATION = "https://api.ncbi.nlm.nih.gov"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov"
_CHR = re.compile(r"^(?:[1-9]|1[0-9]|2[0-2]|X|Y)$", re.I)
# Versioned GRCh38 / GRCh37 primary chromosome accessions, NCBI RefSeq.
_CHROM_VERSIONS = {
    "GRCh38": [11,12,12,12,10,12,14,11,12,11,10,12,11,9,10,10,11,10,10,11,9,11],
    "GRCh37": [10,11,11,11,9,11,13,10,11,10,9,11,10,8,9,9,10,9,9,10,8,10],
}
_XY = {"GRCh38":{"X":11,"Y":10},"GRCh37":{"X":10,"Y":9}}

class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError("NCBI reference service unexpectedly redirected.")

_opener = build_opener(HTTPSHandler(), _NoRedirect())


def chromosome_accession(chromosome, assembly):
    ch = str(chromosome).removeprefix("chr").upper()
    if assembly not in _CHROM_VERSIONS or not _CHR.fullmatch(ch):
        raise ValueError("Use a GRCh38/GRCh37 primary chromosome 1–22, X or Y.")
    if ch in ("X","Y"):
        return f"NC_{23 if ch == 'X' else 24:06d}.{_XY[assembly][ch]}"
    number = int(ch)
    return f"NC_{number:06d}.{_CHROM_VERSIONS[assembly][number-1]}"


def genome_snv(chromosome, genomic_position, transcript_ref, transcript_alt, strand, assembly):
    """Orient selected transcript SNV into genomic-plus allele coordinates."""
    if strand not in (-1, 1) or assembly not in _CHROM_VERSIONS:
        raise ValueError("Invalid transcript strand or reference assembly.")
    if isinstance(genomic_position,bool) or not isinstance(genomic_position,int) or genomic_position < 1:
        raise ValueError("Genomic variant position must be a positive integer.")
    ref,alt = str(transcript_ref).upper(),str(transcript_alt).upper()
    if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT" or ref==alt:
        raise ValueError("Only one-base A/C/G/T substitutions can be queried.")
    if strand == -1:
        ref,alt = str(Seq(ref).complement()),str(Seq(alt).complement())
    accession = chromosome_accession(chromosome,assembly)
    ch=str(chromosome).removeprefix("chr").upper()
    return {
        "assembly":assembly,"chromosome":ch,"genomic_position_1based":genomic_position,
        "genomic_ref":ref,"genomic_alt":alt,"chromosome_refseq":accession,
        "genomic_hgvs":f"{accession}:g.{genomic_position}{ref}>{alt}",
        "spdi":f"{accession}:{genomic_position-1}:{ref}:{alt}",
        "gnomad_query":f"{ch}:{genomic_position}:{ref}:{alt}({assembly})",
    }


def _fetch_json(host, path, params=None):
    if host not in (NCBI_VARIATION,NCBI_EUTILS) or not path.startswith("/") or ".." in path:
        raise ValueError("Unsafe NCBI query destination.")
    url=host+path
    if params: url+="?"+urlencode(params)
    req=Request(url,headers={"Accept":"application/json","User-Agent":"Genotype-PCR-RFLP-Research/0.6"})
    try:
        with _opener.open(req,timeout=13) as response:
            raw=response.read(1_500_001)
            if len(raw)>1_500_000:raise ValueError("Reference response exceeds safe size.")
            return json.loads(raw.decode("utf-8"))
    except HTTPError as exc:
        if exc.code==404:
            return {"_not_found":True}
        raise ValueError(f"NCBI returned HTTP {exc.code} at {path}.") from exc
    except (TimeoutError,URLError,OSError) as exc:
        raise ValueError("NCBI endpoint timed out or could not be reached.") from exc
    except (ValueError,UnicodeError,json.JSONDecodeError) as exc:
        raise ValueError("NCBI returned an invalid or oversized response.") from exc


def query_dbsnp(variant):
    """Allele-based SPDI rsID lookup using the exact assembly chromosome."""
    spdi=variant["spdi"]
    response=_fetch_json(NCBI_VARIATION,
        "/variation/v0/spdi/"+quote(spdi,safe=":")+"/rsids")
    if response.get("_not_found"):
        return {"status":"NOT_FOUND_IN_QUERY","rsids":[],"source":"NCBI Variation Services",
                "note":"No rsID returned for this genomic allele; not evidence of novelty."}
    raw=response.get("rsids")
    if not isinstance(raw,list):
        return {"status":"UNVERIFIED_RESPONSE","rsids":[],"source":"NCBI Variation Services",
                "note":"Unexpected API result format; rsID matching could not be verified."}
    rsids=sorted({"rs"+str(x).removeprefix("rs") for x in raw
                  if isinstance(x,(int,str)) and str(x).removeprefix("rs").isdigit()})
    return {"status":"ALLELE_ASSOCIATED_RSID" if rsids else "NOT_FOUND_IN_QUERY",
            "rsids":rsids[:30],"source":"NCBI Variation Services",
            "note":"An rsID is an identifier for a reference allele record, not a pathogenicity classification."
                   if rsids else "No rsID returned; does not establish a novel variant."}


def query_clinvar(variant):
    """Exact-looking genomic allele query; records remain *candidates*."""
    term=variant["gnomad_query"]
    response=_fetch_json(NCBI_EUTILS,"/entrez/eutils/esearch.fcgi",{
        "db":"clinvar","term":term,"retmode":"json","retmax":10,"tool":"genotype_research_pipeline",
        "email":""  # no user address transmitted
    })
    if response.get("_not_found"):
        return {"status":"UNAVAILABLE","variation_ids":[],"note":"ClinVar query unavailable."}
    result=response.get("esearchresult")
    if not isinstance(result,dict) or not isinstance(result.get("idlist"),list):
        return {"status":"UNVERIFIED_RESPONSE","variation_ids":[],"note":"Invalid ClinVar ESearch response."}
    ids=[str(x) for x in result["idlist"] if str(x).isdigit()]
    total=int(result.get("count",len(ids)))
    return {
        "status":"CANDIDATE_RECORDS_REQUIRE_ALLELE_CHECK" if total else "NOT_FOUND_IN_QUERY",
        "reported_count":total,"variation_ids":ids[:10],
        "links":[f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{i}/" for i in ids[:10]],
        "source":"ClinVar E-utilities",
        "note":("ClinVar search returned candidate records. Verify the EXACT REF/ALT, genome build, "
                "and interpretation on each record before claiming a known variant."
                if total else "No match returned for this exact genomic query; not evidence of novelty."),
    }


def inspect_variant(chromosome, position, ref, alt, strand, assembly):
    """Independent ClinVar/dbSNP checks; one failure does not suppress the other."""
    variant=genome_snv(chromosome,position,ref,alt,strand,assembly)
    results={"variant":variant,"checked_at_utc":datetime.now(timezone.utc).isoformat(),
             "status":"DATABASE_QUERY_ONLY_NOT_CLINICAL_VALIDATION"}
    for label,func in (("dbsnp",query_dbsnp),("clinvar",query_clinvar)):
        try:
            results[label]=func(variant)
        except (ValueError,KeyError,TypeError) as exc:
            results[label]={"status":"QUERY_FAILED","note":str(exc)}
    if results["dbsnp"]["status"]=="ALLELE_ASSOCIATED_RSID":
        results["novelty_assessment"]="Previously catalogued in dbSNP for the queried allele; ClinVar status is separate."
    else:
        results["novelty_assessment"]=(
            "Not established. A negative or failed database lookup NEVER proves that a "
            "variant is novel. Additional sources and literature review are required."
        )
    return results
