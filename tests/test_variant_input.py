"""Tests use artificial DNA only."""
import pytest
from src.variant_input import parse_marked_reference, compare_wild_type_mutant

SYNTHETIC = ("ACGT" * 25) + "[GAG]" + ("TGCA" * 25)

def test_marked_snv():
    wt, mutation, meta = parse_marked_reference(SYNTHETIC, "AAG")
    assert len(wt) == 203
    assert mutation == "G101A"
    assert meta["mutant_full_sequence"][100:103] == "AAG"

def test_paired_sequence():
    wt, mutation, meta = parse_marked_reference(SYNTHETIC, "AAG")
    seq, alt, info = compare_wild_type_mutant(wt, meta["mutant_full_sequence"])
    assert seq == wt and alt == mutation

@pytest.mark.parametrize("motif", ["GAG","AAA","AG","ATGT"])
def test_invalid_change(motif):
    with pytest.raises(ValueError):
        parse_marked_reference(SYNTHETIC, motif)

def test_prose_rejected():
    with pytest.raises(ValueError):
        parse_marked_reference("Subject: sample\nAAA[GAG]CCC","AAG")
