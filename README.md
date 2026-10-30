# Genotype — PCR-RFLP assay designer

Interactive, research-use, in-silico assay design for a single-base substitution in a supplied DNA sequence.

## Features
- Strict single-record FASTA and reference-allele validation
- Primer3 flanking primer design, with clearly labelled provisional fallback
- Differential restriction enzyme screening of the actual PCR amplicon
- Wild-type, heterozygote, and mutant fragment predictions
- Corrected illustrative virtual gel
- Streamlit web app with FASTA upload, interactive tables, PNG, JSON/CSV/ZIP downloads
- CLI single-variant and batch CSV workflow
- Automated unit tests and GitHub Actions CI

## Install and run locally

Requires Python 3.10–3.12 (recommended).

~~~bash
cd /Volumes/Extreme_SSD/Genotype_Pipeline
source venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py
~~~

This local path is on the user's Mac; updated files must first be checked out from the feature branch.

## CLI

~~~bash
python cli.py analyze --seq examples/ttn_exon344_HMERF.fasta --mutation C101T --name TTN_HMERF
python cli.py batch --csv examples/batch_example.csv
python -m pytest -q
~~~

## Coordinate and clinical safety notes

Variant notation `C101T` or `101C>T` means the **101st position within the supplied sequence**. A prefix such as `c.` does not implement HGVS reference transcript mapping; input must be mapped to the FASTA explicitly. A reference-base mismatch produces an error.

A `DESIGN_CANDIDATE` means only an in-silico differential-digest candidate. This software **does not** confirm genome-wide primer uniqueness, amplification, reaction conditions, methylation-sensitive enzyme behavior, experimental genotyping, full HGVS mapping, external restriction-enzyme availability, or validated dCAPS primers. The mismatch-primer fallback is deliberately disabled pending rigorous allele-level implementation. Small restriction fragments may be poorly resolved in agarose gels.

The old Markdown/Word/HTML report modules remain in the repo for backward compatibility, but this feature branch CLI saves JSON and a PNG gel; the Streamlit app exports ZIP/CSV/JSON/PNG. Fuller PDF/Word reporting is a later task.
