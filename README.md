# Genotype Pipeline

Reusable local tool for designing PCR-RFLP assays for any DNA sequence + point mutation.

Given a sequence and a mutation, it automatically:
- Finds commercial restriction enzymes that cut differently in wild-type vs mutant
- Designs flanking primers (Primer3 + mismatch fallback)
- Predicts fragment sizes
- Generates a virtual gel image
- Writes Markdown, Word (.docx), and HTML reports
- Supports batch processing of multiple mutations

Designed to be reusable for any gene / any point mutation.

## Requirements

- Python 3.9 or higher
- Git

**Install all Python packages:**

```bash
pip install biopython primer3-py click pyyaml python-docx matplotlib pandas

Or:
Bashpip install -r requirements.txtInstallationmacOS / LinuxBashgit clone https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline.git
cd Genotype_Any_Point_Mutation_Pipeline
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

Windows (Command Prompt)cmdgit clone https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline.git
cd Genotype_Any_Point_Mutation_Pipeline
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt

Windows (PowerShell)PowerShellgit clone https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline.git
cd Genotype_Any_Point_Mutation_Pipeline
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txtUsageSingle mutationBashpython cli.py analyze \
  --seq examples/ttn_exon344_HMERF.fasta \
  --mutation C101T \
  --name "TTN_HMERF_C30729Y"Batch modeBashpython cli.py batch --csv examples/batch_example.csv

TTN Examples1. HMERF (c.92185C>T, p.Cys30729Tyr)Primers (Pair 1):
Forward: AGGAAAAGTAGTTGGGCCAA (Tm=56.9°C)
Reverse: GTGTAACTTTGGGCTTCACA (Tm=56.2°C)
Product: 180 bp
Virtual Gel:2. TMD / Udd Myopathy (c.107884A>G, p.Glu35962Lys)Primers (Pair 1):
Forward: TCCTCCATGCAAGAGTCCTT (Tm=58.3°C)
Reverse: ACTTTGCCTTCATCAATGCT (Tm=55.9°C)
Product: 170 bp
Virtual Gel:

OutputThe pipeline generates:
Virtual gel image (.png)
Markdown report (.md)
Word report (.docx)
HTML report (.html)
All files are saved in the output/ folder.


AuthorMudasir Rashid
Postdoctoral Associate, Howard University College of Medicine
Google Scholar
ORCID
Loop
text
