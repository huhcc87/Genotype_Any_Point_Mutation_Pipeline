"""Single source of truth for CLI, Streamlit and batch PCR-RFLP analysis."""
from datetime import datetime, timezone
from .sequence import load_sequence, parse_mutation, apply_mutation
from .restriction import find_differential_enzymes
from .primers import design_primers, simple_end_primers, check_primer_quality

DISCLAIMER = (
    "Research-use in-silico prediction only. Validate genome-wide primer specificity, "
    "reference transcript/genome coordinates, enzyme availability and methylation "
    "sensitivity, reaction conditions, and experimental genotype bands."
)


def analyze_assay(seq_input: str, mutation: str, name="Variant", max_enzymes=10):
    sequence = load_sequence(seq_input)
    pos, ref, alt = parse_mutation(mutation, sequence)
    wt, mut = apply_mutation(sequence, pos, ref, alt)
    primers = design_primers(sequence, pos)
    warnings = [DISCLAIMER, "Coordinates are 1-based within the supplied sequence, NOT automatically HGVS or genomic coordinates."]
    if not primers:
        primers = simple_end_primers(sequence, mutation_pos=pos)
        if primers:
            warnings.append("Primer3 found no pair: provisional end primers require Tm, specificity and PCR validation.")
    if not primers:
        warnings.append("No suitable primer pair found; no PCR-RFLP assay can be recommended.")
    candidates = []
    for pair in primers:
        left, right = pair["forward_start"], pair["reverse_start"]
        if not (0 <= left < pos < right < len(sequence)):
            continue
        pair["warnings"] = (
            check_primer_quality(pair["forward"], pair["forward_tm"]) +
            check_primer_quality(pair["reverse"], pair["reverse_tm"])
        )
        amp_wt = wt[left:right + 1]
        amp_mut = mut[left:right + 1]
        enzymes = find_differential_enzymes(amp_wt, amp_mut, max_results=max_enzymes)
        pair["amplicon_start_1based"] = left + 1
        pair["amplicon_end_1based"] = right + 1
        pair["mutation_in_amplicon_1based"] = pos - left + 1
        pair["enzyme_count"] = len(enzymes)
        if enzymes:
            candidates.append((enzymes[0]["resolvable"], enzymes[0]["score"], -len(pair["warnings"]), pair, enzymes))
    if not candidates:
        status = "NO_DIFFERENTIAL_ENZYME" if primers else "NO_PRIMERS"
        warnings.append("No discriminatory restriction digest confirmed for the candidate amplicons.")
        chosen = primers[0] if primers else None
        enzymes = []
    else:
        candidates.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
        chosen, enzymes = candidates[0][3], candidates[0][4]
        status = "DESIGN_CANDIDATE" if enzymes[0]["resolvable"] else "LOW_RESOLUTION_CANDIDATE"
        if not enzymes[0]["resolvable"]:
            warnings.append("Differential fragments may be too small to resolve on standard agarose.")
    return {
        "name": str(name),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "position": pos + 1,
        "ref": ref,
        "alt": alt,
        "mutation": f"{ref}{pos + 1}{alt}",
        "seq_length": len(sequence),
        "coordinate_system": "1-based position within supplied sequence",
        "status": status,
        "primers": primers,
        "selected_primer": chosen,
        "enzymes": enzymes,
        "warnings": warnings,
        "wt_amplicon": wt[chosen["forward_start"]:chosen["reverse_start"] + 1] if chosen else None,
        "mut_amplicon": mut[chosen["forward_start"]:chosen["reverse_start"] + 1] if chosen else None,
    }
