"""Strict local-sequence and single-nucleotide substitution handling."""
from io import StringIO
from pathlib import Path
import re
from Bio import SeqIO

MAX_LENGTH = 2_000_000
DNA = frozenset("ACGT")


def load_sequence(seq_input: str) -> str:
    """Accept a literal DNA sequence, FASTA text, or a local FASTA file path."""
    if not isinstance(seq_input, str) or not seq_input.strip():
        raise ValueError("Provide a DNA sequence or FASTA file.")
    source = seq_input.strip()
    if len(source) < 4096 and "\n" not in source:
        try:
            path = Path(source)
            if path.is_file():
                if path.suffix.lower() not in (".fa", ".fasta", ".fna", ".txt"):
                    raise ValueError("Sequence files must be FASTA (.fa/.fasta/.fna/.txt).")
                source = path.read_text(encoding="utf-8")
        except OSError as exc:
            raise ValueError("Cannot read sequence path.") from exc
    if source.startswith(">"):
        records = list(SeqIO.parse(StringIO(source), "fasta"))
        if len(records) != 1:
            raise ValueError("FASTA input must contain exactly one record.")
        source = str(records[0].seq)
    sequence = re.sub(r"\s+", "", source).upper()
    if not sequence or len(sequence) > MAX_LENGTH:
        raise ValueError(f"Sequence length must be between 1 and {MAX_LENGTH:,} bp.")
    invalid = set(sequence) - DNA
    if invalid:
        raise ValueError("Sequence contains non-ACGT characters: " + "".join(sorted(invalid)))
    return sequence


def parse_mutation(mutation: str, sequence: str = None):
    """Parse a substitution as a 1-based position relative to the supplied sequence.

    Supported: G123A, 123G>A, c.G123A, c.123G>A.
    c. is only input syntax here; true HGVS transcript mapping is NOT inferred.
    Returns (0-based index, reference base, alternate base).
    """
    if not isinstance(mutation, str):
        raise ValueError("Mutation must be a string.")
    value = mutation.strip().upper()
    if value.startswith("C."):
        value = value[2:]
    match = re.fullmatch(r"([ACGT])([1-9]\d*)([ACGT])", value)
    if match:
        ref, position, alt = match.groups()
    else:
        match = re.fullmatch(r"([1-9]\d*)([ACGT])>([ACGT])", value)
        if not match:
            raise ValueError("Use G123A or 123G>A, with position relative to supplied sequence.")
        position, ref, alt = match.groups()
    if ref == alt:
        raise ValueError("Reference and alternate alleles must differ.")
    index = int(position) - 1
    if sequence is not None:
        if index >= len(sequence):
            raise ValueError("Mutation position is outside the supplied sequence.")
        if sequence[index].upper() != ref:
            raise ValueError(
                f"Reference mismatch at position {index + 1}: "
                f"input says {ref}, sequence contains {sequence[index].upper()}."
            )
    return index, ref, alt


def apply_mutation(sequence: str, position: int, ref: str, alt: str):
    """Return WT and single-nucleotide mutant sequences (position is 0-based)."""
    if not isinstance(position, int) or position < 0 or position >= len(sequence):
        raise ValueError("Mutation position is outside the supplied sequence.")
    ref, alt = ref.upper(), alt.upper()
    if ref not in DNA or alt not in DNA or ref == alt:
        raise ValueError("Mutation must substitute one A/C/G/T nucleotide.")
    if sequence[position].upper() != ref:
        raise ValueError(
            f"Reference mismatch at position {position + 1}: "
            f"expected {ref}, found {sequence[position].upper()}."
        )
    wt = sequence.upper()
    return wt, wt[:position] + alt + wt[position + 1:]
