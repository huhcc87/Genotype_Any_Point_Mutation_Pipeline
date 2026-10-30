# Security policy

This is a locally run research prototype for **in-silico PCR-RFLP assay design**, not a clinical or protected-data application.

## Report a vulnerability privately

Use the repository's GitHub **Security → Report a vulnerability** option if available, or privately contact the repository maintainer through an established institutional channel. Please do not submit exploits, unpublished DNA sequences, credentials, participant identifiers, or protected health information in public issues.

## Deployment and data handling

Run Streamlit bound to `127.0.0.1` on a trusted workstation. A public Streamlit deployment is **not** appropriate for patient information without independent institutional security review. Uploaded data are processed in memory and may be present in session state and generated downloads while a session is active. Exported JSON and FASTA contain the sequence, so handle and store outputs accordingly.

Supported uploads are limited to FASTA, TXT, DOCX, and text-based PDF (10 MB). This is an input-size control, **not a full sandbox** for parsing hostile files. Do not accept untrusted DOCX/PDF files from unknown sources until resource-limited parsing isolation is in place.

Dependencies are checked in CI. Future security work should include dependency locks, content limits for decompressed document parsing, limits on server memory/time, and secret scanning.

## Responsible use

Do not provide patient genotype interpretations or clinical decisions. The current application predicts candidate laboratory assays rather than calling a genotype from experimental data.
