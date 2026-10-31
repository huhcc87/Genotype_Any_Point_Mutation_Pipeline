from src.gene_structure import transcript_structure, gene_locus_summary
import pytest

def test_plus_strand_exons_introns_utrs():
    t={"id":"NM_000001.1","strand":1,"seq_region_name":"2",
       "Translation":{"start":103,"end":209},
       "Exon":[{"start":100,"end":105},{"start":200,"end":212}]}
    result=transcript_structure(t)
    s=result["summary"]
    assert s["total_exons"]==2 and s["total_introns"]==1
    assert s["transcript_genomic_span_bp"]==113
    assert s["spliced_transcript_length_bp"]==19
    assert result["introns"][0]["intron_length_bp"]==94
    assert result["junctions"][0]["donor_intronic_positions_5to3"]==[106,107]
    assert result["exons"][0]["utr_5prime_bp"]==3
    assert result["exons"][1]["utr_3prime_bp"]==3

def test_minus_strand_exons_and_junction_orientation():
    t={"id":"NM_000002.1","strand":-1,"seq_region_name":"2",
       "Translation":{"start":103,"end":209},
       "Exon":[{"start":100,"end":105},{"start":200,"end":212}]}
    result=transcript_structure(t)
    assert result["exons"][0]["genomic_start_1based"]==200
    assert result["junctions"][0]["donor_intronic_positions_5to3"]==[199,198]
    assert result["junctions"][0]["acceptor_intronic_positions_5to3"]==[107,106]
    assert result["exons"][0]["utr_5prime_bp"]==3
    assert result["exons"][1]["utr_3prime_bp"]==3

def test_coding_triplets_exclude_terminal_stop_from_protein():
    t={"id":"NM_000003.1","strand":1,"seq_region_name":"2",
       "Translation":{"start":100,"end":108},
       "Exon":[{"start":100,"end":108}]}
    s=transcript_structure(t,"ATGAAATAA")["summary"]
    assert s["coding_triplets_including_stop"]==3
    assert s["predicted_protein_length_aa"]==2
    assert s["terminal_stop_codon"]=="TAA"

def test_gene_locus_is_separate_from_selected_transcript():
    x=gene_locus_summary({"gene_symbol":"GENE1","chromosome":"chr2",
                          "assembly":"GRCh38","gene_locus":{"start_1based":100,"end_1based":200}})
    assert x["gene_span_bp"]==101

def test_invalid_overlapping_exons_rejected():
    t={"strand":1,"Exon":[{"start":1,"end":20},{"start":18,"end":25}]}
    with pytest.raises(ValueError, match="overlap"):
        transcript_structure(t)
