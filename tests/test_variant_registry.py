"""Allele orientation and conservative record evidence tests (no network)."""
from src import variant_registry as v

def test_forward_reverse_genomic_variant():
    plus=v.genome_snv("chr2",100,"C","T",1,"GRCh38")
    assert plus["spdi"]=="NC_000002.12:99:C:T"
    minus=v.genome_snv("chr2",100,"G","A",-1,"GRCh38")
    assert minus["spdi"]=="NC_000002.12:99:C:T"
    assert minus["genomic_hgvs"]=="NC_000002.12:g.100C>T"

def test_no_match_is_not_called_novel(monkeypatch):
    monkeypatch.setattr(v,"_fetch_json",lambda host,path,params=None: {"rsids":[]} if "rsids" in path else {"esearchresult":{"count":"0","idlist":[]}})
    result=v.inspect_variant("chr2",123,"C","T",1,"GRCh38")
    assert result["dbsnp"]["status"]=="NOT_FOUND_IN_QUERY"
    assert result["clinvar"]["status"]=="NOT_FOUND_IN_QUERY"
    assert "not" in result["novelty_assessment"].lower()

def test_clinvar_results_are_only_candidate_records(monkeypatch):
    monkeypatch.setattr(v,"_fetch_json",lambda host,path,params=None: {"esearchresult":{"count":"1","idlist":["123456"]}})
    result=v.query_clinvar(v.genome_snv("chr2",123,"C","T",1,"GRCh38"))
    assert result["status"]=="CANDIDATE_RECORDS_REQUIRE_ALLELE_CHECK"

def test_registry_failure_fails_soft(monkeypatch):
    monkeypatch.setattr(v,"_fetch_json",lambda *a,**kw: (_ for _ in ()).throw(ValueError("offline")))
    result=v.inspect_variant("chr2",123,"C","T",1,"GRCh38")
    assert result["dbsnp"]["status"]=="QUERY_FAILED"
    assert result["clinvar"]["status"]=="QUERY_FAILED"
    assert "Not established" in result["novelty_assessment"]

def test_assembly_difference():
    assert v.chromosome_accession("chr2","GRCh37")=="NC_000002.11"
    assert v.chromosome_accession("chr2","GRCh38")=="NC_000002.12"
