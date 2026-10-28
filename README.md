# Genotype Pipeline

Reusable local tool for designing PCR-RFLP assays for any DNA sequence + point mutation.

Given a sequence and a mutation, it automatically:
- Finds commercial restriction enzymes that cut differently in wild-type vs mutant
- Designs flanking primers
- Predicts fragment sizes
- Generates a virtual gel image
- Writes a Markdown report

## Installation

```bash
git clone https://github.com/huhcc87/Genotype-Any_Point-Mutation-_Pipeline.git
cd Genotype-Any_Point-Mutation-_Pipeline
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## TTN Examples

### 1. HMERF (c.92185C>T, p.Cys30729Tyr)

**Primers (Pair 1):**
- Forward: `AGGAAAAGTAGTTGGGCCAA` (Tm=56.9°C)
- Reverse: `GTGTAACTTTGGGCTTCACA` (Tm=56.2°C)
- Product: 180 bp

**Amplicon size:** 180 bp  
**Expected bands (SetI):**
- Wild-type (undigested): 180 bp
- Wild-type (digested): ~100 bp + ~80 bp
- Mutant (undigested / not cut): 180 bp
- Heterozygote: 180 bp + 100 bp + 80 bp

The mutation destroys the SetI site, so the mutant allele remains uncut (single 180 bp band), while the wild-type allele is cut into two smaller fragments.
**Virtual Gel:**

![TTN HMERF Virtual Gel](docs/TTN_HMERF_gel.png)

### 2. TMD / Udd Myopathy (c.107884A>G, p.Glu35962Lys)

**Primers (Pair 1):**
- Forward: `TCCTCCATGCAAGAGTCCTT` (Tm=58.3°C)
- Reverse: `ACTTTGCCTTCATCAATGCT` (Tm=55.9°C)
- Product: 170 bp
**Amplicon size:** 170 bp  
**Expected bands (AbaSI):**
- Wild-type (undigested): 170 bp
- Wild-type (digested): smaller fragments
- Mutant: differential pattern depending on site loss/gain

AbaSI cuts differently between wild-type and mutant, producing a distinguishable fragment pattern on the gel.
**Virtual Gel:**

![TTN TMD Virtual Gel](docs/TTN_TMD_gel.png)

### Summary Table

| Variant | Mutation | Enzyme | Amplicon | WT Digested | Mutant | Heterozygote |
|---------|----------|--------|----------|-------------|--------|--------------|
| HMERF   | c.92185C>T | SetI | 180 bp | ~100 + 80 bp | 180 bp | 180 + 100 + 80 bp |
| TMD     | c.107884A>G | AbaSI | 170 bp | Cut | Uncut / different | Both patterns |
## Author

**Mudasir Rashid**

- Google Scholar: https://scholar.google.com/citations?user=EjdVSIcAAAAJ&hl=en
- ORCID: https://orcid.org/0000-0001-9342-6243
