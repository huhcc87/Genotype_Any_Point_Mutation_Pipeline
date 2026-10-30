"""Synthetic-only tests. Never include unpublished research DNA in public tests."""
import pytest
from src.variant_input import parse_marked_reference, compare_wild_type_mutant

SYNTHETIC = ("ACGT" * 25) + "[GAG]" + ("TGCA" * 25)

def test_annotated_one_base_change():
    wt, variant, details = parse_marked_reference(SYNTHETIC, "AAG")
    assert len(wt) == 203
    assert variant == "G101A"
    assert wt[100:103] == "GAG"
    assert details["mutant_full_sequence"][100:103] == "AAG"

def test_compare_full_wt_and_mutant():
    wt, variant, details = parse_marked_reference(SYNTHETIC, "AAG")
    seq, derived, meta = compare_wild_type_mutant(wt, details["mutant_full_sequence"])
    assert seq == wt
    assert derived == variant
    assert meta["local_substitution_1based"] == 101

@pytest.mark.parametrize("replacement", ["GAG", "AAA", "AG", "ATGT"])
def test_reject_zero_multiple_or_length_changes(replacement):
    with pytest.raises(ValueError):
        parse_marked_reference(SYNTHETIC, replacement)

def test_reject_prose_and_multiple_marks():
    with pytest.raises(ValueError):
        parse_marked_reference("Subject: sample\nAAA[GAG]CCC", "AAG")
    with pytest.raises(ValueError):
        parse_marked_reference("AAA[GAG]CCC[GAG]TTT", "AAG")

def test_reject_indel_and_multisite_mutant():
    with pytest.raises(ValueError):
        compare_wild_type_mutant("AACCGG", "AACCG")
    with pytest.raises(ValueError):
        compare_wild_type_mutant("AACCGG", "AGTCGG")
