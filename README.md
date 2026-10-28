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

**Python packages** (installed automatically):
- biopython
- primer3-py
- click
- pyyaml
- python-docx
- matplotlib
- pandas

## Installation

### macOS / Linux

```bash
git clone https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline.git
cd Genotype_Any_Point_Mutation_Pipeline
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
