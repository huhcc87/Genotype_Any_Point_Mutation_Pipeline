from src.variant_registry import query_dbsnp, curated_public_reference, genome_snv

def test_documented_dbsnp_response_shape(monkeypatch):
    from src import variant_registry as v
    monkeypatch.setattr(v,"_fetch_json",lambda *args,**kwargs:{"data":{"rsids":[11568350]}})
    record=v.query_dbsnp(v.genome_snv("chr2",189565370,"C","A",1,"GRCh38"))
    assert record["rsids"]==["rs11568350"]

def test_q248h_reference_only_exact_genomic_allele():
    hit=curated_public_reference(genome_snv("chr2",189565370,"C","A",1,"GRCh38"))
    assert hit["coding_hgvs"]=="NM_014585.6:c.744G>T"
    assert hit["rsid"]=="rs11568350"
    assert curated_public_reference(genome_snv("chr2",189565370,"C","T",1,"GRCh38")) is None
