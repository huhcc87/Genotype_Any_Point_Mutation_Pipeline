# Genotype — interactive PCR-RFLP assay designer

A research-use **in-silico PCR-RFLP design** application. It accepts DNA sequences from FASTA, TXT, DOCX and text-based PDFs, designs primer candidates, ranks differential restriction enzymes on the actual primer amplicon, simulates complete digestion and exports multi-format reports.

## Launch locally on macOS

```bash
cd /Volumes/Extreme_SSD/Genotype_Pipeline
git fetch origin
git switch feature/streamlit-validated-pcr-rflp
source venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q
streamlit run app.py
```

The current venv may use Python 3.9; the CI matrix targets 3.10, 3.11, and 3.12. Prefer a new Python 3.11 virtual environment for reproducibility.

## Five-step workflow

1. Select **Try example** for a ready-to-run 200 bp TTN HMERF test sequence (C101T in the local fragment), or upload FASTA/TXT/Word DOCX/text-selectable PDF, or paste a DNA sequence.
2. Enter a substitution in 1-based **local sequence** coordinates, e.g. `C101T`. The tool **does not automatically map HGVS c. coordinates** to genomic/exon sequence.
3. Review selected forward/reverse primer sequences, melting temperatures, colored binding regions, variant location and product size.
4. Compare **design grades** A–D and heuristic scores for potential restriction enzymes. These prioritize visible discriminatory fragments and simple cut patterns; scores are **not** measurements of enzyme specificity, cost, star activity or performance.
5. Select a DNA ladder, visualize virtual WT/heterozygous/mutant digestion lanes, and download **all results** as one ZIP.

## Report bundle

`assay_report.pdf`, `assay_report.docx`, `assay_report.html`, `sequence_highlighted.html`, `virtual_gel.png`, `primer_map.png`, `all_primers.csv`, `ranked_enzymes.csv`, `expected_fragments.csv`, `analysis.json`, `input_reference.fasta`, `README_RESULTS.txt`.

## Supported upload conventions

- FASTA: one or multiple records with standard `>` headers; choose the intended record from the list.
- TXT/DOCX/PDF: provide DNA sequence on separate contiguous lines (ideally at least 60 bases). Multiple blocks can be individually selected.
- Scanned PDFs cannot be read without OCR. Binary Word `.doc` is not supported; save as DOCX.
- Upload limit 10 MB. Only A/C/G/T sequence characters are accepted. Do not upload patient identifiers to publicly hosted deployments.
- Check the **reference allele** at the variant position before analysis. Coordinate system is always local, not automatically transcript HGVS.

## Available marker options

Thermo Fisher GeneRuler 100 bp SM0241 and NEB Quick-Load Purple 100 bp N0551 (manufacturer-derived nominal sizes), plus **generic illustrative** 50 bp and 1 kb presets. Simulated migration is not calibrated and cannot predict precise band brightness, resolution, or digestion yield. For a real gel, use manufacturer instructions for the precise catalogue product.

## CLI

```bash
python cli.py analyze --seq examples/ttn_exon344_HMERF.fasta --mutation C101T --name TTN_HMERF
python cli.py batch --csv examples/batch_example.csv
python -m pytest -q
```

## Scientific limitations

Research use only; not a clinical diagnostic or sample genotype caller. In-silico candidates require independently validated primer uniqueness, lab PCR, fragment visualization, enzyme activity and complete digestion assumptions. Full HGVS/reference accession mapping, experimentally validated dCAPS mismatch primers, genome-wide primer specificity and experimental reaction optimisation are **not implemented**. Do not order reagents based on grading alone.

## Development

App: `app.py`; analysis: `src/assay.py`; input: `src/document_input.py`; ladder metadata: `src/ladders.py`; ranking: `src/ranking.py`; figures: `src/visualization.py`; reports: `src/exports.py`; test: `tests/`. Verify GitHub Actions before merging the feature branch.
