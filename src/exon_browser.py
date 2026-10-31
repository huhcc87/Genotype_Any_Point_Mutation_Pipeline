"""Presentation helpers for transcript-wide coding-codon and exon browsing.

All codon positions are globally 1-based in the selected transcript's coding
sequence. A splice-junction codon may belong to more than one exon.
"""
from collections import Counter


def summarize_coding_exons(transcript, cds, rows, coding_map_function):
    """Summarize coding base and codon coverage for every coding exon."""
    coords, labels = coding_map_function(transcript, cds)
    if len(coords) != len(cds) or len(labels) != len(cds):
        raise ValueError("Inconsistent coding-coordinate map.")
    base_count = Counter(labels)
    result = []
    for exon in sorted(base_count):
        membership = [row for row in rows if exon in row["exon_numbers"]]
        if not membership:
            continue
        result.append({
            "exon_number": exon,
            "coding_bases": base_count[exon],
            "codons_overlapping_exon": len(membership),
            "first_codon": min(row["codon_number"] for row in membership),
            "last_codon": max(row["codon_number"] for row in membership),
            "junction_codons": sum(bool(row["splice_junction_codon"]) for row in membership),
        })
    if sum(row["coding_bases"] for row in result) != len(cds):
        raise ValueError("Exon summary does not account for all coding bases.")
    return result


def filter_codon_rows(rows, exon=None):
    """None means all coding exons; numeric exon filters junction codons too."""
    if exon is None:
        return list(rows)
    if not isinstance(exon, int) or isinstance(exon, bool) or exon < 1:
        raise ValueError("Exon number must be a positive integer.")
    return [row for row in rows if exon in row["exon_numbers"]]


def lookup_codon(rows, codon_number):
    """Retrieve a globally numbered coding codon regardless of exon filter."""
    if not isinstance(codon_number, int) or isinstance(codon_number, bool):
        raise ValueError("Codon number must be an integer.")
    if not 1 <= codon_number <= len(rows):
        raise ValueError("Codon number outside this transcript's CDS.")
    result = rows[codon_number - 1]
    if result["codon_number"] != codon_number:
        raise ValueError("Unexpected codon indexing; refusing reference prediction.")
    return result
