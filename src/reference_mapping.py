"""Map a user-supplied contiguous genomic reference segment to local SNV coordinates.

Important: this does NOT resolve an accession from NCBI, verify a genome build,
or implement spliced transcript/HGVS c. mapping. User-supplied metadata are
trusted only as *claims*; the reference base in the supplied sequence is checked.
"""
import re
from Bio.Seq import Seq
from .sequence import load_sequence

_GENOMIC = re.compile(
    r"(?:(?P<accession>[A-Z]{1,4}_[0-9]+\.[0-9]+):)?"
    r"g\.(?P<position>[1-9][0-9]*)(?P<ref>[ACGT])>(?P<alt>[ACGT])",
    flags=re.IGNORECASE,
)
_ACCESSION = re.compile(r"^(?:NC|NG|NW|NT)_[0-9]+\.[0-9]+$", flags=re.IGNORECASE)


def map_genomic_snv(sequence, notation, accession, region_start, strand="+", assembly="Unspecified"):
    """Map g. position to 1-based local position within contiguous DNA.

    The *region_start* is the 1-based position of the lower-coordinate end of
    the interval on the accession's genomic plus strand, independent of whether
    the supplied DNA sequence is shown forward (+) or reverse-complemented (-).

    Returns local allele text suitable for the existing assay engine plus
    explicit metadata describing the limitation of this mapping.
    """
    dna = load_sequence(sequence)
    acc = str(accession).strip().upper()
    if not _ACCESSION.fullmatch(acc):
        raise ValueError("Use a versioned genomic accession, for example NC_000001.11.")
    if strand not in ("+", "-"):
        raise ValueError("Strand must be + or -.")
    if isinstance(region_start, bool) or not isinstance(region_start, int) or region_start < 1:
        raise ValueError("Region start must be a positive 1-based genomic coordinate.")
    text = str(notation).strip()
    match = _GENOMIC.fullmatch(text)
    if not match:
        raise ValueError(
            "Use NC_000001.11:g.12345C>T or g.12345C>T. "
            "Spliced transcript c./n. variants are not supported in genomic mapping mode."
        )
    notation_acc = match["accession"]
    if notation_acc and notation_acc.upper() != acc:
        raise ValueError("Variant accession does not match the supplied reference accession.")
    position = int(match["position"])
    end = region_start + len(dna) - 1
    if not region_start <= position <= end:
        raise ValueError(f"Genomic coordinate {position} is outside the declared interval {region_start}-{end}.")
    ref, alt = match["ref"].upper(), match["alt"].upper()
    if ref == alt:
        raise ValueError("Reference and alternate alleles must differ.")
    if strand == "+":
        index = position - region_start
        local_ref, local_alt = ref, alt
    else:
        index = end - position
        local_ref = str(Seq(ref).complement())
        local_alt = str(Seq(alt).complement())
    if dna[index] != local_ref:
        raise ValueError(
            f"Reference-base mismatch: at local position {index + 1}, "
            f"expected {local_ref} from {acc}:g.{position}{ref}>{alt} "
            f"({strand} strand), but supplied DNA contains {dna[index]}."
        )
    local = f"{local_ref}{index + 1}{local_alt}"
    return {
        "local_mutation": local,
        "reference_accession": acc,
        "coordinate_system": "genomic plus-strand g. substitution mapped to supplied contiguous region",
        "region_start_1based": region_start,
        "region_end_1based": end,
        "input_sequence_strand": strand,
        "assembly": str(assembly or "Unspecified"),
        "original_genomic_variant": f"{acc}:g.{position}{ref}>{alt}",
        "mapping_status": "LOCAL_BASE_VERIFIED_ONLY",
        "mapping_warning": (
            "The genomic interval, assembly and accession are user supplied; "
            "the app verifies only the base in the submitted DNA fragment. "
            "No external accession retrieval, alignment or HGVS transcript mapping has been performed."
        ),
    }
