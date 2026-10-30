from pathlib import Path
import zipfile
from src.assay import analyze_assay
from src.ranking import rank_enzymes
from src.exports import build_export_bundle
from src.ladders import LADDERS

def test_report_package_example():
    sequence=Path("examples/ttn_exon344_HMERF.fasta").read_text()
    result=analyze_assay(sequence,"C101T","TTN HMERF demo",10)
    result["enzymes"]=rank_enzymes(result["enzymes"])
    enzyme=result["enzymes"][0] if result["enzymes"] else None
    from src.sequence import load_sequence
    zipped,individual=build_export_bundle(result,load_sequence(sequence),
                                          next(iter(LADDERS)),enzyme)
    with zipfile.ZipFile(__import__("io").BytesIO(zipped)) as archive:
        names=set(archive.namelist())
        assert {"assay_report.pdf","assay_report.docx",
                "sequence_highlighted.html","ranked_enzymes.csv",
                "input_reference.fasta","analysis.json"}.issubset(names)
        assert archive.testzip() is None
