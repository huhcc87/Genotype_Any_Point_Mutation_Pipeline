"""Primer3 secondary-structure metrics plus exact repeats in supplied DNA only.

None of these checks is a genome-wide specificity or wet-lab validation result.
"""
from Bio.Seq import Seq
import primer3


def _tm(result):
    return round(float(result.tm), 1) if result.structure_found else None


def _sites(haystack, needle):
    out = []
    start = 0
    while True:
        index = haystack.find(needle, start)
        if index < 0:
            break
        out.append(index + 1)
        start = index + 1
    return out


def assess_primer_pair(pair, reference_sequence):
    """Return JSON-compatible metrics, quality warnings, and scope string."""
    forward, reverse = pair["forward"].upper(), pair["reverse"].upper()
    source = reference_sequence.upper()
    if not source or any(base not in "ACGT" for base in source + forward + reverse):
        raise ValueError("Primer quality checks require A/C/G/T-only DNA.")
    fwd_rc = str(Seq(forward).reverse_complement())
    rev_rc = str(Seq(reverse).reverse_complement())
    # Count exact binding sequences on either strand of supplied fragment.
    forward_sites = {
        (position, "+") for position in _sites(source, forward)
    } | {(position, "-") for position in _sites(source, fwd_rc)}
    reverse_sites = {
        (position, "+") for position in _sites(source, reverse)
    } | {(position, "-") for position in _sites(source, rev_rc)}
    metrics = {
        "forward_hairpin_tm_C": _tm(primer3.calc_hairpin(forward)),
        "reverse_hairpin_tm_C": _tm(primer3.calc_hairpin(reverse)),
        "forward_homodimer_tm_C": _tm(primer3.calc_homodimer(forward)),
        "reverse_homodimer_tm_C": _tm(primer3.calc_homodimer(reverse)),
        "heterodimer_tm_C": _tm(primer3.calc_heterodimer(forward, reverse)),
        "forward_exact_hits_in_supplied_sequence": len(forward_sites),
        "reverse_exact_hits_in_supplied_sequence": len(reverse_sites),
        "specificity_scope": "Exact full-length matches in supplied DNA fragment only; NOT a genome-wide screen",
    }
    notes = []
    for label in (
        "forward_hairpin_tm_C", "reverse_hairpin_tm_C",
        "forward_homodimer_tm_C", "reverse_homodimer_tm_C", "heterodimer_tm_C"
    ):
        value = metrics[label]
        if value is not None and value >= 45:
            notes.append(f"Potential secondary structure: {label}={value} C at Primer3 default chemistry")
    for side in ("forward", "reverse"):
        count = metrics[f"{side}_exact_hits_in_supplied_sequence"]
        if count != 1:
            notes.append(f"{side.capitalize()} primer has {count} exact strand-aware matches in the supplied fragment (expected one).")
    notes.append(
        "Genome-wide specificity and primer-site polymorphism checks were NOT performed. "
        "Verify against the relevant genome/reference using Primer-BLAST or an equivalent tool."
    )
    metrics["quality_warnings"] = notes
    return metrics
