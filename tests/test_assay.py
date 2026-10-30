import pytest
from Bio.Restriction import EcoRI
from src.sequence import load_sequence,parse_mutation,apply_mutation
from src.restriction import digest_fragments,restriction_cuts
from src.gel import gel_fragment_profiles

def test_fasta():
    assert load_sequence(">example\nACGT\nACGT\n")=="ACGTACGT"
    with pytest.raises(ValueError): load_sequence("ACGTXYZ")

@pytest.mark.parametrize("mutation",["A2G","2A>G","c.A2G","c.2A>G"])
def test_notation(mutation):
    assert parse_mutation(mutation,"CAC")==(1,"A","G")

@pytest.mark.parametrize("mutation",["A0G","A2A","A2GXXX","A99G","G2A"])
def test_reject(mutation):
    with pytest.raises(ValueError): parse_mutation(mutation,"CAC")

def test_bounds():
    with pytest.raises(ValueError): apply_mutation("ACG",-1,"A","G")

def test_fragments():
    assert digest_fragments(100,[25,60])==[40,35,25]
    assert gel_fragment_profiles(100,[25],[])["Heterozygote"]==[75,25,100]

def test_coordinates():
    assert restriction_cuts(EcoRI,"AAAAGAATTCTTT")==[5]
