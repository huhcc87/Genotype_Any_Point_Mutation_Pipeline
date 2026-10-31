from src.exon_browser import summarize_coding_exons, filter_codon_rows, lookup_codon
from src.gene_lookup import coding_map, codon_table
import pytest

def test_all_exons_and_splice_junctions():
    transcript={"strand":1,"Translation":{"start":100,"end":204},
                "Exon":[{"start":100,"end":103},{"start":200,"end":204}]}
    cds="ATGGAATAA"
    rows=codon_table(transcript,cds)
    summary=summarize_coding_exons(transcript,cds,rows,coding_map)
    assert len(rows)==3
    assert sum(x["coding_bases"] for x in summary)==9
    assert len(filter_codon_rows(rows,None))==3
    assert len(filter_codon_rows(rows,1))==2
    assert len(filter_codon_rows(rows,2))==2
    assert lookup_codon(rows,2)["exon_numbers"]==[1,2]

def test_direct_codon_validation():
    rows=[{"codon_number":1,"reference_codon":"ATG","exon_numbers":[1]}]
    assert lookup_codon(rows,1)["reference_codon"]=="ATG"
    with pytest.raises(ValueError):
        lookup_codon(rows,2)
