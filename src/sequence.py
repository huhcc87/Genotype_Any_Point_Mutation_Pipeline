"""Strict DNA input parsing. Web inputs are NEVER interpreted as filesystem paths."""
from io import StringIO
from pathlib import Path
import re
from Bio import SeqIO

MAX_LENGTH = 2_000_000
DNA = frozenset("ACGT")
ALLOWED_FASTA_SUFFIXES = frozenset({".fa", ".fasta", ".fna", ".txt"})


def load_sequence(seq_input: str) -> str:
    """Parse only literal A/C/G/T or a single FASTA record, never a file path.

    This behavior is deliberate: web form data must not grant local file read
    access. The CLI uses the separate explicit read_sequence_file function.
    """
    if not isinstance(seq_input, str) or not seq_input.strip():
        raise ValueError("Provide a DNA sequence or FASTA text.")
    source = seq_input.strip()
    if len(source) > MAX_LENGTH * 2 + 4096:
        raise ValueError("Input text exceeds the supported size limit.")
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


def read_sequence_file(path: str) -> str:
    """Read an explicitly requested local FASTA/TXT file for trusted CLI use."""
    file = Path(path).expanduser()
    if file.suffix.lower() not in ALLOWED_FASTA_SUFFIXES:
        raise ValueError("CLI file must be .fa, .fasta, .fna or .txt.")
    if not file.is_file():
        raise ValueError(f"Not a readable sequence file: {file}")
    if file.stat().st_size > MAX_LENGTH * 2 + 4096:
        raise ValueError("Sequence file exceeds the size limit.")
    return load_sequence(file.read_text(encoding="utf-8-sig"))


def parse_mutation(mutation: str, sequence: str = None):
    """Parse A123G or 123A>G relative to 1-based supplied DNA coordinates.

    A c. prefix is accepted as simplified input notation; HGVS transcript
    mapping is NOT inferred and an accession is not resolved.
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
            raise ValueError("Use G123A or 123G>A in supplied-sequence coordinates.")
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
    """Return WT and mutant DNA for an explicitly validated single substitution."""
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
