import pytest
from src.variant_input import parse_marked_reference, compare_wild_type_mutant

MSH3 = ("TCCTGAAAAGCAGACAGGTTTCTAACACTTTAGAAATAGAGTACATACATACTCCTGAGTGTATCATCTTTACTTTTTCCATGTTTATGCTGTGTTATAGCCATAACTGGGGAAATACATTTTTTCTGTAACATTATATTTGTATTTGTTTTTAGGGAGTGCAGCCTGCCACAGGCGAGGTTGTGTTTGATAGTTTCCAGGACTCTGCTTCTCGTTCA[GAG]CTAGAAACCCGGATGTCAAGCCTGCAGCCAGTAGAGCTGCTGCTTCCTTCGGCCTTGTCCGAGCAAACAGAGGCGCTCATCCACAGAGCCACATCTGTTAG")

def test_msh3_e413k_annotated_example():
    wt,mutation,info=parse_marked_reference(MSH3,"AAG")
    assert len(wt)==322
    assert mutation=="G219A"
    assert wt[218:221]=="GAG"
    assert info["mutant_full_sequence"][218:221]=="AAG"

def test_whole_sequence_comparison():
    wt,mutation,meta=parse_marked_reference(MSH3,"AAG")
    seq,mut,details=compare_wild_type_mutant(wt,meta["mutant_full_sequence"])
    assert seq==wt and mut==mutation

@pytest.mark.parametrize("motif",["GAG","AAA","AG","ATGT"])
def test_reject_other_than_one_snv(motif):
    with pytest.raises(ValueError):
        parse_marked_reference(MSH3,motif)

def test_reject_prose_and_multiple_marks():
    with pytest.raises(ValueError):
        parse_marked_reference("Subject: exon 8\nAAA[GAG]CCC","AAG")
    with pytest.raises(ValueError):
        parse_marked_reference("AAA[GAG]CCC[GAG]TTT","AAG")

def test_paired_sequences_reject_indels_and_multiple_changes():
    with pytest.raises(ValueError):
        compare_wild_type_mutant("AACCGG","AACCG")
    with pytest.raises(ValueError):
        compare_wild_type_mutant("AACCGG","AGTCGG")
