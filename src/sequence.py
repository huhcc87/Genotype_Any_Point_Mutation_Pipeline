from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re

def load_sequence(seq_input: str) -> str:
    """Load sequence from string or FASTA file."""
    seq_input = seq_input.strip()
    if seq_input.endswith((".fa", ".fasta", ".txt")):
        record = SeqIO.read(seq_input, "fasta")
        return str(record.seq).upper()
    # Clean plain sequence
    seq = re.sub(r"[^ATCGatcg]", "", seq_input)
    return seq.upper()

def apply_mutation(sequence: str, position: int, ref: str, alt: str) -> tuple[str, str]:
    """
    Apply a point mutation.
    position is 0-based.
    Returns (wildtype, mutant)
    """
    sequence = sequence.upper()
    ref = ref.upper()
    alt = alt.upper()

    if sequence[position] != ref:
        raise ValueError(
            f"Reference base mismatch at position {position}. "
            f"Expected {ref}, found {sequence[position]}"
        )

    mutant = sequence[:position] + alt + sequence[position+1:]
    return sequence, mutant

def parse_mutation(mutation: str, sequence: str = None):
    """
    Parse mutation strings like:
    - c.G1237A
    - G1237A
    - 1237G>A
    Returns (0-based_position, ref, alt)
    """
    mutation = mutation.strip().upper().replace("C.", "")
    
    # Pattern: 1237G>A or G1237A
    match = re.match(r"([ACGT])?(\d+)([ACGT])>?([ACGT])?", mutation)
    if not match:
        raise ValueError(f"Cannot parse mutation: {mutation}")

    groups = match.groups()
    if groups[0] and groups[2]:  # G1237A style
        ref, pos, alt = groups[0], int(groups[1]) - 1, groups[2]
    elif groups[1] and groups[2] and groups[3]:  # 1237G>A style
        pos, ref, alt = int(groups[1]) - 1, groups[2], groups[3]
    else:
        raise ValueError(f"Cannot parse mutation: {mutation}")

    return pos, ref, alt
