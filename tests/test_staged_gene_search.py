from src import gene_lookup as g

def test_staged_gene_search_avoids_expanded_gene(monkeypatch):
    requests=[]
    def fake(path, params=None, assembly='GRCh38'):
        requests.append((path,params))
        if path.startswith("/lookup/symbol/"):
            return {"object_type":"Gene","species":"homo_sapiens","assembly_name":"GRCh38","id":"ENSG00000155657","display_name":"TTN","canonical_transcript":"ENST00000591111.2"}
        if path.startswith("/overlap/id/"):
            return [{"id":"ENST00000591111","Parent":"ENSG00000155657","biotype":"protein_coding","strand":-1,"seq_region_name":"2"},
                    {"id":"ENST00000000002","Parent":"ENSG00000000000","biotype":"protein_coding","strand":1,"seq_region_name":"2"}]
        raise AssertionError(path)
    monkeypatch.setattr(g,"_get_json",fake)
    found=g.search_gene("TTN")
    assert len(found["transcripts"])==1
    assert found["transcripts"][0]["is_canonical"]
    assert requests[0][1] is None
    assert requests[1][1] == {"feature":"transcript","biotype":"protein_coding"}

def test_transcript_details_are_lazy(monkeypatch):
    def fake(path,params=None):
        assert path=="/lookup/id/ENST00000591111"
        assert params=={"expand":1}
        return {"object_type":"Transcript","id":"ENST00000591111","Parent":"ENSG00000155657","biotype":"protein_coding","Translation":{"start":1,"end":3},"Exon":[{"start":1,"end":3}],"seq_region_name":"2"}
    monkeypatch.setattr(g,"_get_json",fake)
    result=g.get_transcript_details("ENST00000591111","ENSG00000155657")
    assert result["id"]=="ENST00000591111"
