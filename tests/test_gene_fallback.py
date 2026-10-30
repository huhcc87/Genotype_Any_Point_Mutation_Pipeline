"""Reference-discovery fallback tests use mocked API data only."""
from src import gene_lookup as g
import pytest

def test_fallback_when_ensembl_gene_endpoint_fails(monkeypatch):
    def primary(path, params=None, assembly="GRCh38"):
        raise ValueError("Ensembl HTTP 500")
    monkeypatch.setattr(g, "_get_json", primary)
    monkeypatch.setattr(g, "_mygene_lookup", lambda name: {
        "gene_id":"ENSG00000155657", "gene_symbol":name,
        "assembly":"GRCh38","transcripts":[{"id":"ENST00000591111"}],
        "reference_provider":"MyGene.info (discovery only)","lookup_warning":"Not verified",
        "query":name,
    })
    result = g.search_gene("TTN")
    assert result["gene_symbol"] == "TTN"
    assert result["reference_provider"].startswith("MyGene.info")
    assert "HTTP 500" in result["lookup_warning"]

def test_no_fallback_cross_assembly(monkeypatch):
    def primary(path, params=None, assembly="GRCh38"):
        raise ValueError("Ensembl HTTP 500")
    monkeypatch.setattr(g, "_get_json", primary)
    monkeypatch.setattr(g, "_mygene_lookup", lambda name: pytest.fail("GRCh37 must not use GRCh38 fallback"))
    with pytest.raises(ValueError, match="GRCh37"):
        g.search_gene("TTN", assembly="GRCh37")

def test_fallback_fails_closed(monkeypatch):
    monkeypatch.setattr(g, "_get_json", lambda *a, **kw: (_ for _ in ()).throw(ValueError("HTTP 500")))
    monkeypatch.setattr(g, "_mygene_lookup", lambda x: (_ for _ in ()).throw(ValueError("unavailable")))
    with pytest.raises(ValueError, match="fallback also failed"):
        g.search_gene("TTN")
