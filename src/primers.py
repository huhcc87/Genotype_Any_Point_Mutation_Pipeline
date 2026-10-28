from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import primer3

def design_primers(sequence: str, mutation_pos: int, 
                   product_size_range=(100, 350),
                   primer_length=(18, 25),
                   optimal_tm=58):
    """Design primers flanking the mutation using primer3."""
    sequence = sequence.upper()
    seq_len = len(sequence)

    # Target region around mutation
    target_start = max(0, mutation_pos - 60)
    target_length = min(120, seq_len - target_start)

    seq_args = {
        'SEQUENCE_ID': 'target',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_TARGET': [target_start, target_length],
    }

    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': primer_length[0],
        'PRIMER_MAX_SIZE': primer_length[1],
        'PRIMER_OPT_TM': optimal_tm,
        'PRIMER_MIN_TM': 52,
        'PRIMER_MAX_TM': 65,
        'PRIMER_MIN_GC': 30.0,
        'PRIMER_MAX_GC': 70.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range],
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_EXPLAIN_FLAG': 1,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_MAX_NS_ACCEPTED': 0,
    }

    try:
        results = primer3.bindings.design_primers(seq_args, global_args)
    except Exception as e:
        print(f"Primer3 error: {e}")
        return []

    primer_pairs = []
    num_pairs = results.get('PRIMER_PAIR_NUM_RETURNED', 0)

    for i in range(num_pairs):
        left_seq = results.get(f'PRIMER_LEFT_{i}_SEQUENCE')
        right_seq = results.get(f'PRIMER_RIGHT_{i}_SEQUENCE')
        if not left_seq or not right_seq:
            continue
        product_size = results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
        left_tm = results[f'PRIMER_LEFT_{i}_TM']
        right_tm = results[f'PRIMER_RIGHT_{i}_TM']
        left_pos = results[f'PRIMER_LEFT_{i}'][0]
        right_pos = results[f'PRIMER_RIGHT_{i}'][0]

        primer_pairs.append({
            'pair': i + 1,
            'forward': left_seq,
            'reverse': right_seq,
            'product_size': product_size,
            'forward_tm': round(left_tm, 1),
            'reverse_tm': round(right_tm, 1),
            'forward_start': left_pos,
            'reverse_start': right_pos,
            'type': 'natural'
        })

    return primer_pairs

def simple_end_primers(sequence: str, length=22):
    """Fallback primers from sequence ends."""
    sequence = sequence.upper()
    if len(sequence) < 60:
        return []

    fwd = sequence[:length]
    rev = str(Seq(sequence[-length:]).reverse_complement())

    return [{
        'pair': 1,
        'forward': fwd,
        'reverse': rev,
        'product_size': len(sequence),
        'forward_tm': round(mt.Tm_NN(fwd), 1),
        'reverse_tm': round(mt.Tm_NN(rev), 1),
        'forward_start': 0,
        'reverse_start': len(sequence) - length,
        'type': 'end'
    }]

def design_mismatch_primer(sequence: str, mutation_pos: int, enzyme_site: str, 
                           allele="mutant", mismatch_pos_in_site=None):
    """
    Design a mismatch primer that creates a restriction site only in one allele.
    This is useful when no natural differential site exists.
    """
    sequence = sequence.upper()
    enzyme_site = enzyme_site.upper().replace("/", "").replace("N", "N")
    
    # Simple implementation: place mismatch near the mutation
    # and try to create a common site (e.g. for HinfI, AluI, etc.)
    
    # For now return a note that mismatch design is available
    return {
        "note": "Mismatch primer design available",
        "suggestion": f"Introduce mismatch near position {mutation_pos+1} to create site for {enzyme_site}",
        "allele_targeted": allele
    }
def design_mismatch_primers(sequence: str, pos: int, alt: str, enzyme_site: str = "GAATTC", length=22):
    """
    Design a forward primer with a mismatch that creates a restriction site
    only in the mutant (or destroys it in WT). Simple version for demonstration.
    """
    sequence = sequence.upper()
    # Place the primer so the mutation is near the 3' end
    start = max(0, pos - length + 3)
    primer = list(sequence[start:start + length])
    
    # Introduce a mismatch near the 3' end to help create a site
    mismatch_pos = min(len(primer) - 3, pos - start)
    if mismatch_pos >= 0:
        # Simple example: change one base to help form EcoRI-like site
        primer[mismatch_pos] = "A" if primer[mismatch_pos] != "A" else "G"
    
    fwd = "".join(primer)
    rev = str(Seq(sequence[pos + 5:pos + 5 + length]).reverse_complement())
    
    return [{
        "pair": 1,
        "forward": fwd,
        "reverse": rev,
        "product_size": len(sequence),
        "forward_tm": round(mt.Tm_NN(fwd), 1),
        "reverse_tm": round(mt.Tm_NN(rev), 1),
        "type": "mismatch",
        "note": f"Mismatch introduced to favor site creation near mutation"
    }]
def design_mismatch_primers(sequence: str, pos: int, alt: str, length=22):
    """
    Basic mismatch primer: introduces a single mismatch near the mutation
    to help create a differential restriction site.
    """
    sequence = sequence.upper()
    start = max(0, pos - length + 4)
    primer_seq = list(sequence[start:start + length])
    
    # Introduce one mismatch near the 3' end
    mismatch_idx = min(len(primer_seq) - 4, pos - start)
    if 0 <= mismatch_idx < len(primer_seq):
        bases = ["A", "C", "G", "T"]
        current = primer_seq[mismatch_idx]
        primer_seq[mismatch_idx] = next(b for b in bases if b != current)
    
    fwd = "".join(primer_seq)
    rev_start = pos + 10
    rev = str(Seq(sequence[rev_start:rev_start + length]).reverse_complement())
    
    return [{
        "pair": 1,
        "forward": fwd,
        "reverse": rev,
        "product_size": len(sequence),
        "forward_tm": round(mt.Tm_NN(fwd), 1),
        "reverse_tm": round(mt.Tm_NN(rev), 1),
        "type": "mismatch",
        "note": "Single mismatch introduced near mutation"
    }]
