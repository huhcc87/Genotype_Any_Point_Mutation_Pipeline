# Genotype Pipeline

Reusable local tool for designing PCR-RFLP assays for any DNA sequence + point mutation.

Given a sequence and a mutation, it automatically:
- Finds commercial restriction enzymes that cut differently in wild-type vs mutant
- Designs flanking primers
- Predicts fragment sizes
- Generates a virtual gel image
- Writes a Markdown report

Designed to be reusable for any gene / any point mutation.

## Installation

```bash
git clone https://github.com/huhcc87/Genotype-Any_Point-Mutation-_Pipeline.git
cd Genotype-Any_Point-Mutation-_Pipeline
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
