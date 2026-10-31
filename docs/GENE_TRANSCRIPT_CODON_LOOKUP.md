# Search by gene, transcript, exon, and coding codon

This research module is provided on the feature branch `feature/gene-transcript-exon-codon-lookup`. Internet access is required to retrieve *public reference* annotations from Ensembl. The application does not upload patient data to Ensembl.

## Workflow

1. Under **Input source**, select **Search gene / transcript / exon / codon**.
2. Enter a human gene symbol or an Ensembl stable **ENSG** identifier; press **Find gene and transcripts**.
3. Select the protein-coding **ENST** transcript. The app displays accession, version (if returned), coding protein length and exon count. Codon numbering can differ between transcripts.
4. Press **Load CDS and exon/codon map**. Select a coding exon. The table lists codon numbers, the three-base WT codon, amino acid, CDS base positions and genomic coordinates. If a codon crosses an exon junction, it is flagged.
5. Select one codon; then select a **single-base mutant codon** and review the predicted amino-acid effect.
6. Press **Load reference genomic DNA for selected variant**. The app obtains a GRCh38 genomic DNA region around the variant, orients it according to transcript strand and confirms the corresponding WT base.
7. Press **Analyze sequence** to calculate PCR primer candidates, differential restriction patterns, a virtual gel, and the downloadable report package.

## Scope and important caveats

- This implementation is limited to **human protein-coding Ensembl transcripts** on GRCh38 and **one-nucleotide codon substitutions**. Chromosomes 1–22, X, Y and MT are recognized, although nonstandard mitochondrial translation is not implemented, so mitochondrial coding interpretation should **not** be relied upon.
- Annotation is retrieved from the **current Ensembl REST release**. Transcript versions are displayed where available; the method does not pin a release or download a frozen reference dataset.
- The workflow reconstructs coding position from exon/translation coordinates and cross-checks the number of mapped bases against the retrieved CDS. It validates the target genomic base but does **not** verify all genomic exon bases or run a genome-wide primer specificity scan.
- Different transcript isoforms can assign different codon and exon numbers to the same genomic position. Always confirm the intended clinical/research reference.
- This output is **reference-based assay design**, not experimental genotyping of an unknown sample or a diagnostic interpretation.
- Ensembl uses public gene annotations. Uploaded reference sequences are not transmitted by this module.
- Reference API failure, inconsistent coding sequence lengths and a reference-base mismatch must stop analysis rather than silently guess.

### General example

Select a protein-coding transcript and an exon, then choose a codon from that transcript. The app shows the reference codon and possible single-nucleotide substitutions. No particular gene or variant is preselected.

## Source and scientific method

- [Ensembl REST documentation](https://rest.ensembl.org/)
- Ensembl lookup by symbol/ID, sequence by transcript ID and sequence by human genomic region.
- Genome assembly: GRCh38.

## Browse all coding exons and codons

After **Load CDS and exon/codon map**, the app now shows an **Overview of every coding exon** with coding-base counts, codon ranges and junction-codon counts. Choose **All coding exons** to display every reference coding codon, or select an individual exon to filter the table. Use **Go directly to coding codon / amino-acid number** to locate any 1-based codon in the entire selected transcript, independent of the exon filter. A codon spanning two exons appears in both filtered exon views but only once in the all-codons table. These are **coding** exons: untranslated exons or noncoding portions of exons are not represented in the codon view.

## Full gene structure and splice junction report

After selecting a transcript and loading its coding sequence, the application now displays two **separate** annotation scopes:

- **Gene-level location:** chromosome, genome assembly, gene identifier, source-provided start/end coordinates and whole-gene genomic span in bp/kb, when available. A gene locus can contain multiple isoforms and should not be inferred from a single transcript. Cytogenetic band appears only if independently provided, rather than guessed from chromosome coordinates.
- **Transcript-specific structure:** genomic start/end and span, plus/minus strand, all exon and intron counts, spliced transcript length, coding-exon count, coding region length, 5′ and 3′ UTR lengths, total intronic length, coding triplets including any terminal stop codon, and predicted protein amino-acid count excluding a terminal stop codon.

Open **All exons, introns and splice-junction coordinates** to view every exon including UTR regions, every intervening intron, and each theoretical splice donor/acceptor intronic boundary (1-based genomic coordinates). Exon numbers progress in transcript 5′→3′ order, even on a negative-strand gene. Junctions are derived from the exon model; GT/AG motifs have **not** been sequence-verified. Click **Download gene and transcript structure (JSON)** for the complete result.

**Note:** The number of exons in the coding view can differ from the total transcript exon count. In different isoforms, the number and length of exons/introns can change. Predicted coding triplets can include a stop codon that does not encode an amino acid. All values are specific to the selected reference assembly and annotated transcript.
