"""Conservative PCR primer design using Primer3; no unvalidated dCAPS claims."""
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import primer3


def check_primer_quality(primer: str, tm: float) -> list[str]:
    if not primer:
        return ["Empty primer"]
    gc = 100 * (primer.count("G") + primer.count("C")) / len(primer)
    warnings = []
    if not 40 <= gc <= 60:
        warnings.append(f"GC content {gc:.1f}% outside preferred 40–60%")
    if not 55 <= tm <= 65:
        warnings.append(f"Tm {tm:.1f} °C outside preferred 55–65 °C")
    if any(x not in "ACGT" for x in primer):
        warnings.append("Non-ACGT nucleotide present")
    return warnings


def design_primers(sequence: str, mutation_pos: int,
                   product_size_range=(100, 350), primer_length=(18, 25),
                   optimal_tm=60):
    """Design primers flanking a single-base target; coordinates refer to full sequence."""
    sequence = sequence.upper()
    if not 0 <= mutation_pos < len(sequence):
        raise ValueError("Mutation outside sequence.")
    lo = max(0, mutation_pos - 500)
    hi = min(len(sequence), mutation_pos + 501)
    window = sequence[lo:hi]
    target = mutation_pos - lo
    if len(window) < product_size_range[0]:
        return []
    seq_args = {
        "SEQUENCE_ID": "variant_region",
        "SEQUENCE_TEMPLATE": window,
        "SEQUENCE_TARGET": [target, 1],
    }
    globals_ = {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": int(primer_length[0]),
        "PRIMER_MAX_SIZE": int(primer_length[1]),
        "PRIMER_OPT_TM": optimal_tm,
        "PRIMER_MIN_TM": 54,
        "PRIMER_MAX_TM": 66,
        "PRIMER_MIN_GC": 30,
        "PRIMER_MAX_GC": 70,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_size_range)],
        "PRIMER_NUM_RETURN": 8,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY_TH": 47,
        "PRIMER_MAX_SELF_END_TH": 47,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47,
        "PRIMER_PAIR_MAX_COMPL_END_TH": 47,
    }
    try:
        result = primer3.bindings.design_primers(seq_args, globals_)
    except (OSError, ValueError, RuntimeError) as exc:
        raise ValueError(f"Primer3 failed: {exc}") from exc
    pairs = []
    for i in range(result.get("PRIMER_PAIR_NUM_RETURNED", 0)):
        left = result[f"PRIMER_LEFT_{i}"][0] + lo
        right = result[f"PRIMER_RIGHT_{i}"][0] + lo
        size = right - left + 1
        if not left < mutation_pos < right or size != result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]:
            continue
        fwd, rev = result[f"PRIMER_LEFT_{i}_SEQUENCE"], result[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        pairs.append({
            "pair": len(pairs) + 1,
            "forward": fwd,
            "reverse": rev,
            "forward_tm": round(result[f"PRIMER_LEFT_{i}_TM"], 1),
            "reverse_tm": round(result[f"PRIMER_RIGHT_{i}_TM"], 1),
            "product_size": size,
            "forward_start": left,
            "reverse_start": right,
            "type": "primer3",
        })
    return pairs


def simple_end_primers(sequence: str, mutation_pos=None, length=22):
    """Provisional fallback, never guaranteed to amplify specifically."""
    sequence = sequence.upper()
    if len(sequence) < 60 or len(sequence) > 800:
        return []
    if mutation_pos is not None and not (length <= mutation_pos < len(sequence) - length):
        return []
    fwd = sequence[:length]
    rev = str(Seq(sequence[-length:]).reverse_complement())
    return [{
        "pair": 1, "forward": fwd, "reverse": rev,
        "product_size": len(sequence),
        "forward_tm": round(mt.Tm_NN(fwd), 1),
        "reverse_tm": round(mt.Tm_NN(rev), 1),
        "forward_start": 0,
        "reverse_start": len(sequence) - 1,
        "type": "provisional-end",
    }]


def design_mismatch_primers(sequence: str, pos: int, alt: str, **kwargs):
    """Fail closed: the previous implementation produced unverified dCAPS primers.

    Future dCAPS searches must model the primer-incorporated mismatch and
    prove differential restriction cleavage of both resulting PCR products.
    """
    return []
