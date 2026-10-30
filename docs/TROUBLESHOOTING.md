# Troubleshooting — Genotype PCR-RFLP

[One-code installation](INSTALLATION.md) · [User guide](USER_GUIDE.md) · [Repository](../README.md)

## `ModuleNotFoundError: No module named 'src.visualization'`

This happened on an earlier commit because the visualization module was not committed. It is fixed on branch `feature/streamlit-validated-pcr-rflp`. Update to the latest feature branch. If you installed using the one-code local launcher, run that same code again; the script checks for local changes before updating. **Do not delete your existing research folders or use `git reset --hard`.**

## A test fails or the app fails to start

The one-code launcher deliberately stops on test failure. Inspect the last traceback, verify that you are in a folder with `app.py`, `src/` and `requirements.txt`, and check the [GitHub Actions](https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline/actions) result for your commit. The app was validated in CI using Python 3.10, 3.11 and 3.12. The Mac starter uses Python 3.11.

If the existing `.venv` uses another Python version, **rename** it for preservation rather than deleting it, and rerun setup to create a 3.11 environment.

## I cannot open `localhost:8501`

Read the Terminal output and use the specific **Local URL** printed there. If another Streamlit application uses port 8501, the new instance may use **8502**. Keep the Terminal window running. In a local setup, `localhost` only works on the computer running the process.

## `git` or `python3.11` is not found

Install Git using [git-scm.com](https://git-scm.com/downloads/mac). For Python use [Python's macOS downloads](https://www.python.org/downloads/macos/) or, if Homebrew is installed, the launcher installs `python@3.11` for you. Re-run the entire quick-start code.

## Uploaded PDF/Word file is not accepted

Only text-selectable PDF, DOCX, TXT and FASTA are supported; PDFs consisting only of scanned images require external text extraction. Ensure a clear contiguous A/C/G/T-only DNA block of at least 60 bp. For maximum reliability, create a standard FASTA file and upload it instead.

## Reference base mismatch or mutation out of range

The app numbers **within the selected sequence starting at 1**. Verify the REF base, the coordinate and the sequence orientation. HGVS c. and genomic coordinates are **not mapped** into short fragments automatically.

## No primers or restriction enzymes found

That can be a correct outcome. Ensure adequate flanking DNA and the right mutation coordinate, try an alternate validated target region, and remember some substitutions lack a usable natural differential restriction site. Mismatch-primer design remains disabled pending rigorous verification.

## The download button is missing or report creation fails

First run the analysis, then choose a restriction enzyme and ladder, and click **Prepare all reports**. It can produce JSON/CSV for analyses without an appropriate digest. If report generation raises an exception, share a **sanitized** traceback via GitHub Issues and include the current commit ID.

## Research and safety limitations

A successful software test confirms technical behavior, **not** primer specificity or assay performance. Do not use this tool as a clinical genotype caller. Check primers against a genomic reference, confirm reference mapping, restriction buffer, gel percentages, product sizes, enzyme sensitivity to DNA methylation, and the laboratory result.
