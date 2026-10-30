from src import gene_lookup as g
import pytest

def test_coding_exons_on_plus_strand():
    t={"strand":1,"Translation":{"start":104,"end":202},
       "Exon":[{"start":100,"end":106},{"start":200,"end":210}]}
    # CDS spans positions 104,105,106,200,201,202.
    coords, exons=g.coding_map(t,"GAGAAA")
    assert coords==[104,105,106,200,201,202]
    assert exons==[1,1,1,2,2,2]
    rows=g.codon_table(t,"GAGAAA")
    assert rows[0]["reference_codon"]=="GAG"
    assert rows[1]["exon_numbers"]==[2]

def test_splice_junction_codon():
    t={"strand":1,"Translation":{"start":104,"end":200},
       "Exon":[{"start":100,"end":105},{"start":200,"end":210}]}
    rows=g.codon_table(t,"GAG")
    assert rows[0]["splice_junction_codon"] is True
    assert rows[0]["genomic_positions"]==[104,105,200]

def test_negative_strand_coding_map():
    t={"strand":-1,"Translation":{"start":100,"end":202},
       "Exon":[{"start":100,"end":101},{"start":200,"end":203}]}
    coords,exons=g.coding_map(t,"GAGAA")
    assert coords==[202,201,200,101,100]
    assert exons==[1,1,1,2,2]

def test_cds_mapping_mismatch_rejected():
    t={"strand":1,"Translation":{"start":10,"end":12},
       "Exon":[{"start":10,"end":12}]}
    with pytest.raises(ValueError,match="spans"):
        g.coding_map(t,"GGGGGG")

def test_alt_codons():
    options=g.alternative_codons("GAG")
    assert len(options)==9
    assert any(x["mutant_codon"]=="AAG" and x["mutant_amino_acid"]=="K" for x in options)

def test_genomic_allele_and_negative_orientation(monkeypatch):
    # Gene minus-strand CDS GAG, genomic plus strand complementary CTC
    transcript={"id":"ENST00000000001","strand":-1,"seq_region_name":"1",
                "Translation":{"start":300,"end":302},
                "Exon":[{"start":300,"end":302}]}
    def fake(path, params=None):
        assert path.startswith("/sequence/region/human/")
        return {"seq":"A"*300+"CTC"+"A"*300}
    monkeypatch.setattr(g,"_get_json",fake)
    data=g.fetch_genomic_assay(transcript,"GAG",1,"AAG",300)
    assert data["mutation"]=="G301A"
    assert data["sequence"][300]=="G"
    assert data["metadata"]["genomic_position"]==302

def test_gene_lookup_rejects_bad_query():
    with pytest.raises(ValueError):
        g.search_gene("https://example.com/")
