"""Explicit WT/alternate allele input with an unambiguous single-SNV outcome.

Modes:
1. Annotated reference DNA e.g. AC...[GAG]...T, alternate motif AAG.
2. Equal-length full WT and mutant sequences differing at exactly one site.

These derive a *local 1-based* mutation position and do not infer HGVS
transcript coordinates or independently verify external reference accessions.
"""
import re
from .sequence import load_sequence

MAX_MOTIF_LENGTH = 60

def parse_marked_reference(text: str, alternate_motif: str):
    """Return (WT sequence, local variant text, provenance) from [WT] -> ALT."""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Paste a wild-type DNA sequence with exactly one [marked] WT motif.")
    source = text.strip()
    if source.startswith(">"):
        lines = source.splitlines()
        if len(lines) < 2 or any(line.startswith(">") for line in lines[1:]):
            raise ValueError("Annotated FASTA must contain exactly one record.")
        source = "\n".join(lines[1:])
    if source.count("[") != 1 or source.count("]") != 1:
        raise ValueError("Mark exactly one WT region using [GAG], and enter replacement AAG separately.")
    sequence_text = re.sub(r"\s+", "", source).upper()
    match = re.fullmatch(r"([ACGT]*)\[([ACGT]{1,60})\]([ACGT]*)", sequence_text)
    if not match:
        raise ValueError(
            "Only A/C/G/T and one [WT motif] are allowed in the sequence field. "
            "Remove email headings, arrows, labels, and commentary from the sequence box."
        )
    prefix, motif, suffix = match.groups()
    alt = str(alternate_motif).strip().upper()
    if not re.fullmatch(r"[ACGT]{1,60}", alt):
        raise ValueError("Alternate motif must contain only A/C/G/T, e.g. AAG.")
    if len(motif) != len(alt):
        raise ValueError("WT and alternate motifs must have equal lengths; insertions/deletions are unsupported.")
    differences = [i for i, (a, b) in enumerate(zip(motif, alt)) if a != b]
    if len(differences) != 1:
        raise ValueError(
            f"Found {len(differences)} changed bases in the marked motif. "
            "This PCR-RFLP tool currently supports exactly one substituted base."
        )
    local_idx = len(prefix) + differences[0]
    wt = load_sequence(prefix + motif + suffix)
    variant = f"{motif[differences[0]]}{local_idx + 1}{alt[differences[0]]}"
    metadata = {
        "input_mode": "annotated_wild_type_and_alternate_motif",
        "wild_type_motif": motif,
        "mutant_motif": alt,
        "marked_motif_start_1based": len(prefix) + 1,
        "local_substitution_1based": local_idx + 1,
        "variant": variant,
        "variant_coordinate_system": "within supplied WT fragment; not HGVS transcript position",
        "mutant_full_sequence": wt[:local_idx] + alt[differences[0]] + wt[local_idx + 1:],
    }
    return wt, variant, metadata


def compare_wild_type_mutant(wild_type: str, mutant: str):
    """Require identical sequence lengths with precisely one changed nucleotide."""
    wt = load_sequence(wild_type)
    alt = load_sequence(mutant)
    if len(wt) != len(alt):
        raise ValueError("Complete WT and mutant sequences must have equal lengths; indels are unsupported.")
    positions = [i for i, (x, y) in enumerate(zip(wt, alt)) if x != y]
    if len(positions) != 1:
        raise ValueError(
            f"WT and mutant sequences differ at {len(positions)} positions. "
            "Please submit sequences with exactly one nucleotide substitution."
        )
    index = positions[0]
    variant = f"{wt[index]}{index+1}{alt[index]}"
    return wt, variant, {
        "input_mode": "full_wild_type_and_full_mutant_sequences",
        "local_substitution_1based": index+1,
        "variant": variant,
        "variant_coordinate_system": "within supplied WT fragment; not HGVS transcript position",
        "mutant_full_sequence": alt,
    }
