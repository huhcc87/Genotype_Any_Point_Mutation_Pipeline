# Complete user guide — Genotype PCR-RFLP

[One-code local installation](INSTALLATION.md) · [Troubleshooting](TROUBLESHOOTING.md) · [Repository home](../README.md)

## What this tool does

This is a research-use **candidate assay designer** for single-base substitutions. It accepts a reference DNA sequence and a local substitution position, uses Primer3 to propose PCR primers, searches for restriction enzymes with different predicted wild-type/mutant fragment patterns **within the selected PCR amplicon**, ranks candidate restriction enzymes, draws a schematic virtual gel, and exports results.

It does **not** genotype unknown samples, validate enzyme activity or perform transcript/genome coordinate mapping.

## Step 1 — Launch and choose an input

Install and start Streamlit from GitHub using the [single Mac code block](INSTALLATION.md). Open the local URL displayed in Terminal.

| Input method | Accepted types | How to prepare |
| --- | --- | --- |
| **Try example** | Built-in FASTA | Choose Try example and click Analyze example. |
| **Upload file** | `.fa`, `.fasta`, `.fna` | Standard FASTA header `>identifier` with DNA underneath. Multiple records are selectable. |
| **Upload file** | `.txt` | Include uninterrupted DNA on its own line(s), rather than mixed with descriptive prose. |
| **Upload file** | `.docx` | Put the DNA in a paragraph or table, as consecutive A/C/G/T bases. |
| **Upload file** | `.pdf` | Use a *text-selectable* PDF containing a readable DNA block. |
| **Paste sequence** | Plain DNA or single-record FASTA | Paste into the text field. |

The document uploader accepts up to **10 MB** and looks for DNA sequences at least **60 bp** long. It does not support image-only scanned PDFs, older binary Word `.doc` files, automatic OCR, compressed FASTQ/VCF input, or arbitrary DNA files containing ambiguous IUPAC bases. Verify every extracted base before analysis.

## Alternative input: wild-type and mutant alleles

If you already know the **wild-type codon and mutant codon**, you do not need to calculate the nucleotide position manually.

1. Under **Input source**, choose **WT marked motif → mutant motif**.
2. Paste *only* the reference DNA, marking one codon with square brackets; for example `ACCT[GAG]CTG`. Paste the entire target fragment with flanks, not only the short codon.
3. Enter the mutant motif `AAG`.
4. The app checks for exactly **one** nucleotide substitution and calculates the local 1-based variant position.
5. Click **Analyze sequence**, then inspect primers, restriction patterns, and downloadable reports.

Alternatively select **Complete WT and mutant sequences**, paste both entire sequences of equal length, and the app locates the single difference. It rejects insertions/deletions or multiple substitutions.

For a synthetic illustration, use `ACGT[GAG]TTAC` and the mutant motif `AAG`. In a real assay, supply sufficiently long DNA flanks and verify its reference independently.

**Do not include email titles, gene descriptions, arrows, or prose in the DNA entry box.** Use the assay-name field for labels.

## Step 2 — Set the assay and variant

Set an **Assay name**, such as `Gene_Study_Variant`. Enter a substitution, e.g. **C101T** or **101C>T**, where position 101 refers to the **101st base of the selected sequence**. This is *not necessarily* position 101 in a reference transcript or genome. Prefix `c.` is accepted as text notation only and **does not implement HGVS coordinate mapping**. Enter the REF base precisely as it appears in the uploaded sequence. If it is inconsistent, analysis stops with an error.

The maximum restriction enzyme slider limits the number displayed; it does not increase scientific validation.

## Step 3 — Try the one-click example

On first launch, choose **Try example** in the sidebar. The repository's `examples/ttn_exon344_HMERF.fasta` provides a short reference fragment, and the app begins with local substitution **C101T**. Click **Analyze example**.

After the analysis you can inspect primer candidates, a highlighted DNA position diagram, ranked differential enzymes, and gel simulations. The example is **illustrative**: it is not an independently validated clinical TTN assay.

## Step 4 — Understand each result

| Result tab | Contents | Interpretation |
| --- | --- | --- |
| **Primers & highlighted DNA** | 5′→3′ forward and reverse sequences, melting temperatures, amplicon length, colored primer-binding regions, highlighted WT base | Primers are *candidates*; off-target specificity and real amplification are unverified. |
| **Restriction enzyme grades** | Enzyme, site, predicted WT/mutant fragment sizes, heuristic grade A–D and 0–100 design score | The first choice is a computational priority, **not** a verified best vendor enzyme. |
| **Virtual gel & ladder** | Predicted WT, heterozygote and mutant digests, optional DNA marker sets, expected fragment sizes **in bp** | Only simulated complete digests; gel migration, relative band intensity and resolving power vary experimentally. |
| **Download complete reports** | Button to **Prepare all reports** and complete ZIP / individual outputs | Select enzyme and ladder first, then prepare and download reports. |
| **Limitations** | Scientific limitations and assay warnings | Required checks before purchasing or using an assay. |

**WT** means *wild-type reference allele*; **bp** means *base pairs*, not a molecular-weight unit. A heterozygote carries one WT and one alternate allele. Actual band separations and enzyme cleavage may differ due to gel conditions, methylation, incomplete digestion and fragment sizes below resolution.

Available ladders: **Thermo Scientific GeneRuler 100 bp (SM0241)**, **NEB Quick-Load Purple 100 bp (N0551)**, **generic 50 bp illustration**, and **generic 1 kb illustration**. Vendor-specific presets reflect nominal marker sizes; other presets are not manufacturer product claims. Verify current catalog and buffer instructions with manufacturers.

## Step 5 — Download outputs

Click **Prepare all reports**, then **Download ALL results (ZIP)**. The ZIP includes:

| File | What it contains |
| --- | --- |
| `assay_report.pdf` | Printable assay summary and figures |
| `assay_report.docx` | Editable Word report |
| `assay_report.html` | Browser-readable summary; extract ZIP first to keep linked images |
| `sequence_highlighted.html` | Base-by-base reference with green forward, blue reverse and yellow variant highlights |
| `virtual_gel.png` | Schematic genotype lanes with selected DNA ladder |
| `primer_map.png` | Position diagram of primer regions and substitution |
| `all_primers.csv` | Primer sequences, Tm and coordinates |
| `ranked_enzymes.csv` | Differential enzyme design scores and predicted fragments |
| `expected_fragments.csv` | Individual predicted fragments for WT, heterozygote and mutant |
| `analysis.json` | Structured data, parameters and selected enzyme/ladder |
| `input_reference.fasta` | Input reference sequence used for analysis |
| `README_RESULTS.txt` | Instructions and limitations for exported data |

If no suitable primer or discriminating restriction enzyme exists, a **negative result is valid**; the program should not claim an assay has been validated.

## Data protection

Use the local app for research data where feasible. **Do not upload identifiable patient information** to public deployments or public GitHub issues. This software is **not** a clinically validated genotyping system.

## Troubleshooting

See [TROUBLESHOOTING.md](TROUBLESHOOTING.md), including module errors, Python versions, FASTA/PDF parsing and the browser URL. For a reproducible issue, provide Python version, Git commit, exact error, and **non-sensitive synthetic test data**.
