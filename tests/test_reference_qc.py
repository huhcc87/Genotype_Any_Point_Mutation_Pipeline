import pytest
from src.reference_mapping import map_genomic_snv
from src.primer_qc import assess_primer_pair

def test_forward_coordinate_mapping():
    result = map_genomic_snv("AACGT", "NC_000001.11:g.1001A>G",
                              "NC_000001.11", 1000, "+", "GRCh38")
    assert result["local_mutation"] == "A2G"
    assert result["mapping_status"] == "LOCAL_BASE_VERIFIED_ONLY"

def test_reverse_strand_coordinate_mapping():
    # Genomic plus strand A at position 102 corresponds to T on reverse-complement input.
    result = map_genomic_snv("TTG", "g.102A>G", "NC_000001.11", 100, "-", "GRCh38")
    assert result["local_mutation"] == "T1C"

def test_mismatch_rejected():
    with pytest.raises(ValueError, match="mismatch"):
        map_genomic_snv("AACGT", "g.1001C>T", "NC_000001.11", 1000)

def test_transcript_coordinates_not_misrepresented():
    with pytest.raises(ValueError, match="Spliced transcript"):
        map_genomic_snv("AACGT", "NM_000001.1:c.2A>G", "NC_000001.11", 1000)

def test_unsupported_accession_rejected():
    with pytest.raises(ValueError, match="versioned"):
        map_genomic_snv("AACGT", "g.1001A>G", "NM_001", 1000)

def test_primer_metrics_are_local_not_genomewide():
    pair={"forward":"AACCGGTTAACCGGTTAACC","reverse":"TTGGCCAATTGGCCAATTGG"}
    out=assess_primer_pair(pair,"AACCGGTTAACCGGTTAACC"+"G"*45+"CCAATTGGCCAATTGGCCAA")
    assert out["specificity_scope"].startswith("Exact full-length")
    assert "heterodimer_tm_C" in out
    assert any("Genome-wide specificity" in x for x in out["quality_warnings"])
