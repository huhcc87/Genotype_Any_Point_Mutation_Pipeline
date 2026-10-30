"""Regression checks against untrusted local-path and stale-input hazards."""
from pathlib import Path
import pytest
from src.sequence import load_sequence, read_sequence_file
from src.document_input import extract_sequences

def test_untrusted_text_never_reads_local_path(tmp_path):
    secret = tmp_path / "private.fa"
    secret.write_text(">secret\n" + "ACGT" * 50)
    with pytest.raises(ValueError):
        load_sequence(str(secret))
    assert read_sequence_file(str(secret)) == "ACGT" * 50

def test_document_upload_size_limit():
    with pytest.raises(ValueError, match="10 MB"):
        extract_sequences("input.txt", b"A" * (10 * 1024 * 1024 + 1))

def test_document_extension_rejected():
    with pytest.raises(ValueError):
        extract_sequences("input.py", b"A" * 200)

def test_fasta_input_still_supported():
    assert load_sequence(">example\nACGT\n") == "ACGT"
