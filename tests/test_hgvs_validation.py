import pytest
from src.hgvs_validation import parse_coding_snv,validate_coding_snv,validate_batch

TX={"id":"NM_000000.1","strand":-1,"seq_region_name":"2","Translation":{"start":100,"end":108},"Exon":[{"start":100,"end":108}]}
CDS="ATGGAATAA"

def test_legacy_and_standard_snv():
    assert parse_coding_snv("c.4G>A")==parse_coding_snv("c.G4A")==(4,"G","A")

def test_strand_and_consequence():
    result=validate_coding_snv(TX,CDS,"c.4G>A","p.E2K")
    assert result["reference_codon"]=="GAA"
    assert result["alternate_codon"]=="AAA"
    assert result["genomic_position"]==105
    assert result["genomic_ref"]=="C"
    assert result["genomic_alt"]=="T"

def test_mismatched_protein_fails():
    with pytest.raises(ValueError,match="Protein mismatch"):
        validate_coding_snv(TX,CDS,"c.4G>A","p.E5K")

def test_invalid_reference_base_fails():
    with pytest.raises(ValueError,match="Reference mismatch"):
        validate_coding_snv(TX,CDS,"c.4C>A")

def test_batch_each_failure_is_individual():
    records=validate_batch(TX,CDS,"c.4G>A,p.E2K\nc.5C>A,p.E2K")
    assert records[0]["status"]=="REFERENCE_AND_PROTEIN_VERIFIED"
    assert records[1]["status"]=="MISMATCH_OR_UNSUPPORTED"
