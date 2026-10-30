"""Offline fixtures for UCSC RefSeq coordinates and source separation."""
import pytest
from src import ucsc_reference as u
from src.gene_lookup import coding_map

def test_refseq_coordinates_plus():
    row={"name":"NM_000001.2","name2":"GENE1","chrom":"chr2","strand":"+",
         "exonCount":2,"exonStarts":"100,200,","exonEnds":"106,206,",
         "txStart":100,"txEnd":206,"cdsStart":101,"cdsEnd":205}
    result=u._normalize_refseq(row,"GENE1","chr2")
    assert result["id"]=="NM_000001.2"
    assert result["annotation"]["Translation"]=={"start":102,"end":205}
    assert result["annotation"]["Exon"][0]=={"start":101,"end":106}

def test_refseq_reject_mismatched_gene():
    row={"name":"NM_000001.2","name2":"GENE1","chrom":"chr2","strand":"+",
         "exonCount":1,"exonStarts":"100,","exonEnds":"106,",
         "txStart":100,"txEnd":106,"cdsStart":100,"cdsEnd":106}
    assert u._normalize_refseq(row,"OTHER","chr2") is None

def test_refseq_cds_assembled_from_exons(monkeypatch):
    annotation={"source":"UCSC RefSeq Curated","strand":1,
        "refseq_track":{"chrom":"chr2","txStart":100,"txEnd":210,"cdsStart":100,"cdsEnd":210},
        "Exon":[{"start":101,"end":103},{"start":208,"end":210}]}
    dna="ATG"+("C"*104)+"TAA"
    monkeypatch.setattr(u,"fetch_refseq_region",lambda *a,**k:dna)
    assert u.get_refseq_cds(annotation)=="ATGTAA"
    coords,exons=coding_map({**annotation,"Translation":{"start":101,"end":210}}, "ATGTAA")
    assert coords==[101,102,103,208,209,210]

def test_refseq_assembly_not_mixed():
    with pytest.raises(ValueError):
        u.fetch_refseq_region("chr2",100,200,assembly="GRCh36")
