"""Guided PCR-RFLP design UI: document upload, ranked enzymes, gel, complete reports."""
from pathlib import Path
import tempfile
import hashlib
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from src.assay import analyze_assay
from src.restriction import find_differential_enzymes
from src.sequence import apply_mutation
from src.document_input import extract_sequences
from src.sequence import load_sequence
from src.reference_mapping import map_genomic_snv
from src.variant_input import parse_marked_reference, compare_wild_type_mutant
from src.gene_lookup import search_gene, get_transcript_details, get_cds, codon_table, alternative_codons, fetch_genomic_assay
from src.ucsc_reference import search_refseq_gene, get_refseq_cds, fetch_refseq_assay
from src.exon_browser import summarize_coding_exons, filter_codon_rows, lookup_codon
from src.gene_structure import gene_locus_summary, transcript_structure
from src.variant_registry import inspect_variant
from src.hgvs_validation import validate_batch, validate_coding_snv
from src.ranking import rank_enzymes
from src.ladders import LADDERS,ladder_notes,recommended_ladders
from src.visualization import render_ladder_gel,render_sequence_figure,sequence_markup
from src.exports import build_export_bundle

st.set_page_config(page_title="Genotype | PCR-RFLP",page_icon="🧬",layout="wide")
SAMPLE=Path(__file__).resolve().parent/"examples"/"ttn_exon344_HMERF.fasta"
st.title("🧬 Genotype | PCR-RFLP assay designer")
st.caption("Sequence → primers → ranked restriction enzymes → DNA gel → downloadable reports")
st.warning("Research use only: in-silico predictions are not validated laboratory assays or clinical genotype calls.")

with st.expander("User instructions — how to upload and interpret an assay",expanded=True):
    st.markdown("""
**1. Choose an input.** Use the built-in TTN example, paste a DNA sequence, or upload FASTA,
TXT, Microsoft Word DOCX, or selectable-text PDF. DOCX and PDF uploads are searched for
DNA blocks; choose the correct extracted sequence. Scanned PDFs require text extraction/OCR first.

**2. Enter the mutation.** Use C101T or 101C>T for position 101 of the exact uploaded sequence.
All positions are **1-based local reference positions**. HGVS transcript/genomic coordinates
are not automatically converted. Verify the WT base and the reference accession.

**3. Inspect the output.** Primer sequences, predicted amplicon, restriction-enzyme rankings,
WT and mutant digestion patterns, and a virtual gel are shown. An A grade denotes
computational design priority, NOT experimentally verified enzyme efficiency.

**4. Choose a ladder.** WT means wild-type (reference) and bp means base pairs.
Select Thermo or NEB manufacturer-based marker presets, or an illustrative generic marker.

**5. Download.** Export a ZIP with PDF, DOCX, HTML, figures, CSV, JSON,
highlighted DNA sequence and input FASTA. Test in the laboratory before use.
""")
with st.sidebar:
    st.header("Assay setup")
    mode=st.radio("Input source",["Try example","Search gene / transcript / exon / codon","Upload file","Paste sequence","WT marked motif → mutant motif","Complete WT and mutant sequences"])
    name=st.text_input("Assay name","Variant")
    mutation=st.text_input("Mutation (local sequence position)","C101T")
    max_enzymes=st.slider("Maximum restriction enzymes",3,50,20)
    variant_mode=st.radio("Variant coordinates",["Local fragment (1-based)","Genomic g. in declared interval"])
    reference_accession=""
    genomic_start=1
    genomic_strand="+"
    genome_assembly="Unspecified"
    if variant_mode.startswith("Genomic"):
        reference_accession=st.text_input("Versioned genomic accession","NC_000001.11")
        genomic_start=st.number_input("Genomic interval start (1-based)",min_value=1,value=1,step=1)
        genomic_strand=st.selectbox("Submitted sequence strand",["+","-"])
        genome_assembly=st.text_input("Genome assembly (metadata only)","GRCh38")
        st.caption("Enter g.123C>T or NC_000001.11:g.123C>T above. The interval and accession are not independently fetched.")
    st.divider()
    st.caption("Accepted: .fa, .fasta, .fna, .txt, .docx, .pdf (10 MB maximum).")
    st.caption("For gene search, WT marked, and paired-sequence modes, the mutation is derived from the selected reference/codon; the sidebar Mutation field is ignored.")
    st.caption("Do not upload identifiable patient data to any publicly hosted application.")

sequence=None
variant_metadata=None
input_error=None
if mode=="Search gene / transcript / exon / codon":
    st.subheader("Search human gene → transcript → exon → codon")
    st.caption("Ensembl REST reference lookup. Internet required; no mixing of reference assemblies.")
    lookup_assembly=st.selectbox("Reference assembly",["GRCh38","GRCh37"])
    lookup_provider=st.selectbox("Annotation source",["UCSC RefSeq Curated (alternative)","Ensembl REST"],help="UCSC uses versioned NM_ RefSeq transcript IDs; Ensembl uses ENST IDs. Both use the selected genome assembly.")
    identifier=st.text_input("Human gene symbol or Ensembl gene ID","")
    if st.button("Find gene and transcripts"):
        try:
            with st.spinner("Retrieving gene annotation from Ensembl..."):
                found=search_refseq_gene(identifier,assembly=lookup_assembly) if lookup_provider.startswith("UCSC") else search_gene(identifier,assembly=lookup_assembly)
            st.session_state["lookup_gene"]=found
            st.session_state.pop("lookup_cds",None)
            st.session_state.pop("lookup_assay",None)
        except ValueError as exc:
            st.error(str(exc))
    found=st.session_state.get("lookup_gene")
    if found and (found.get("query","").upper()!=identifier.strip().upper() or found.get("assembly")!=lookup_assembly or (lookup_provider.startswith("UCSC") != found.get("reference_provider","").startswith("UCSC"))):
        found=None
        st.session_state.pop("lookup_gene",None)
        st.session_state.pop("lookup_cds",None)
        st.session_state.pop("lookup_assay",None)
    if found:
        st.success(f"{found['gene_symbol']} | {found['gene_id']} | {found['assembly']}")
        st.caption("Gene discovery source: "+found.get("reference_provider","Unknown"))
        if found.get("lookup_warning"):
            st.warning(found["lookup_warning"])
        st.markdown("#### Gene-level location")
        gene_info=gene_locus_summary(found)
        st.dataframe(pd.DataFrame([gene_info]),use_container_width=True,hide_index=True)
        st.caption("Gene span describes the annotated gene locus, not the length of every transcript. Cytoband is shown only if supplied by a verified source.")
        txs=found["transcripts"]
        tx_index=st.selectbox("Choose protein-coding transcript",list(range(len(txs))),
            format_func=lambda i: f"{txs[i]['id']}.{txs[i]['version'] or '?'} | {txs[i]['label']} | exon/codon details on selection" + (" | canonical" if txs[i]["is_canonical"] else ""),
            key="lookup_transcript")
        tx=txs[tx_index]
        if st.button("Load CDS and exon/codon map"):
            try:
                with st.spinner("Mapping coding bases onto genomic exon coordinates..."):
                    annotation=tx["annotation"] if lookup_provider.startswith("UCSC") else get_transcript_details(tx["id"],found["gene_id"],assembly=lookup_assembly)
                    cds=get_refseq_cds(annotation,assembly=lookup_assembly) if lookup_provider.startswith("UCSC") else get_cds(tx["id"],assembly=lookup_assembly)
                    rows=codon_table(annotation,cds)
                st.session_state["lookup_cds"]={"transcript_id":tx["id"],"cds":cds,"rows":rows,"annotation":annotation}
                st.session_state.pop("lookup_assay",None)
            except ValueError as exc:
                st.error(str(exc))
        loaded=st.session_state.get("lookup_cds")
        if loaded and loaded["transcript_id"]==tx["id"]:
            rows=loaded["rows"]
            structure=transcript_structure(loaded["annotation"],loaded["cds"],rows)
            overview=structure["summary"]
            st.markdown("#### Complete selected-transcript structure")
            st.dataframe(pd.DataFrame([overview]),use_container_width=True,hide_index=True)
            st.caption("All exons include UTRs; coding exons are a separate subset. Coordinates are 1-based inclusive.")
            with st.expander("All exons, introns and splice-junction coordinates",expanded=False):
                st.markdown("**All transcript exons (including untranslated regions)**")
                st.dataframe(pd.DataFrame(structure["exons"]),use_container_width=True,hide_index=True)
                st.markdown("**All transcript introns**")
                st.dataframe(pd.DataFrame(structure["introns"]),use_container_width=True,hide_index=True)
                st.markdown("**Splice donor and acceptor boundaries**")
                st.dataframe(pd.DataFrame(structure["junctions"]),use_container_width=True,hide_index=True)
                st.caption("Donor/acceptor positions are predicted from exon structure; splice motifs have not been sequence-verified.")
            st.download_button("Download gene and transcript structure (JSON)",
                __import__("json").dumps({"gene":gene_info,"structure":structure},indent=2),
                file_name="gene_structure_annotation.json",mime="application/json")
            st.caption(f"Validated CDS: {len(loaded['cds'])} bases, {len(rows)} codons. Exon numbering follows the selected transcript.")
            with st.expander("Find or batch-check coding HGVS variants (c. and p.)", expanded=False):
                st.caption("Use the selected transcript accession. Example format: c.1237G>A,p.E413K. Unmatched coordinates and amino-acid changes are explicitly flagged.")
                single_hgvs=st.text_input("Coding DNA variant (HGVS)",placeholder="c.1237G>A")
                single_protein=st.text_input("Expected protein change (optional)",placeholder="p.E413K")
                if single_hgvs.strip():
                    try:
                        verified=validate_coding_snv(loaded["annotation"],loaded["cds"],single_hgvs,single_protein or None)
                        st.success(f"{verified['coding_hgvs']} | {verified['protein_hgvs']} | codon {verified['codon_number']} | exon(s) {verified['exon_numbers']} | chr{verified['chromosome']}:{verified['genomic_position']}{verified['genomic_ref']}>{verified['genomic_alt']}")
                        if st.button("Use verified codon for assay"):
                            st.session_state["lookup_codon_number"]=int(verified["codon_number"])
                            st.session_state["lookup_target_mutant_codon"]=verified["alternate_codon"]
                            st.rerun()
                    except ValueError as exc:
                        st.error(str(exc))
                batch_text=st.text_area("Batch verification: one c. variant per line, optionally followed by comma and p. change",height=110,placeholder="c.1237G>A,p.E413K\\nc.1397G>A,p.S466N")
                if st.button("Verify all submitted variants"):
                    batch=validate_batch(loaded["annotation"],loaded["cds"],batch_text)
                    st.session_state["hgvs_batch"]=batch
                if st.session_state.get("hgvs_batch"):
                    batch=st.session_state["hgvs_batch"]
                    st.dataframe(pd.DataFrame(batch),hide_index=True,use_container_width=True)
                    st.download_button("Download HGVS reference-check table (CSV)",pd.DataFrame(batch).to_csv(index=False),file_name="hgvs_reference_checks.csv")
            exon_summary=summarize_coding_exons(loaded["annotation"],loaded["cds"],rows,__import__("src.gene_lookup",fromlist=["coding_map"]).coding_map)
            with st.expander("Overview of every coding exon",expanded=True):
                st.dataframe(pd.DataFrame(exon_summary),use_container_width=True,hide_index=True)
                st.caption("A codon can cross an exon boundary. Its codon number is unchanged across the whole transcript.")
            exon_options=["All coding exons"]+[str(entry["exon_number"]) for entry in exon_summary]
            exon_choice=st.selectbox("Show codons from",exon_options,key="lookup_exon_filter")
            chosen_exon=None if exon_choice=="All coding exons" else int(exon_choice)
            exon_rows=filter_codon_rows(rows,chosen_exon)
            st.caption(f"Showing {len(exon_rows):,} of {len(rows):,} coding codons. A junction codon may appear in two individual exon views.")
            if st.session_state.get("lookup_row_scope") != (tx["id"], exon_choice):
                st.session_state["lookup_row_scope"]=(tx["id"], exon_choice)
            displayed=pd.DataFrame(exon_rows)[["codon_number","reference_codon","reference_amino_acid","exon_numbers","cds_positions_1based","genomic_positions","splice_junction_codon"]].reset_index(drop=True)
            selected_event=st.dataframe(displayed,use_container_width=True,hide_index=True,height=380,
                                        selection_mode="single-row",on_select="rerun",key=f"lookup_codon_table_{tx['id']}_{exon_choice}")
            selected_rows=selected_event.selection.rows
            if selected_rows:
                row_number=int(exon_rows[selected_rows[0]]["codon_number"])
                if st.session_state.get("lookup_last_selected_row")!=(tx["id"],exon_choice,row_number):
                    st.session_state["lookup_codon_number"]=row_number
                    st.session_state["lookup_last_selected_row"]=(tx["id"],exon_choice,row_number)
                    st.rerun()
            codon_number=int(st.number_input("Select coding codon / amino-acid number (or click a table row)",
                min_value=1,max_value=len(rows),value=1,step=1,key="lookup_codon_number"))
            codon_row=lookup_codon(rows,codon_number)
            st.info(f"Selected codon {codon_number}: {codon_row['reference_codon']} → amino acid {codon_row['reference_amino_acid']}, exon(s) {', '.join(map(str,codon_row['exon_numbers']))}.")
            if chosen_exon is not None and chosen_exon not in codon_row["exon_numbers"]:
                st.caption("Selected codon is outside this exon filter. Choose All coding exons to display its row.")
            variants=alternative_codons(codon_row["reference_codon"])
            alt_index=st.selectbox("Choose one-base mutant codon",list(range(len(variants))),
                format_func=lambda i: f"{codon_row['reference_codon']} → {variants[i]['mutant_codon']} | {codon_row['reference_amino_acid']}{codon_number}{variants[i]['mutant_amino_acid']}",
                key=f"lookup_alt_{tx['id']}_{codon_number}")
            target_alt=st.session_state.pop("lookup_target_mutant_codon",None)
            if target_alt:
                idx=next((i for i,v in enumerate(variants) if v["mutant_codon"]==target_alt),None)
                if idx is not None:
                    st.session_state[f"lookup_alt_{tx['id']}_{codon_number}"]=idx
                    st.rerun()
            alternative=variants[alt_index]
            st.success(f"Selected reference-based change: {codon_row['reference_codon']} → {alternative['mutant_codon']} | {codon_row['reference_amino_acid']}{codon_number}{alternative['mutant_amino_acid']}")
            st.caption("For example, ATG → CTG | M1L means methionine at amino acid 1 is predicted to change to leucine. This is a selected hypothetical substitution, not a sample genotype.")
            registry_variant_key=(lookup_assembly,tx["id"],codon_number,alternative["mutant_codon"])
            if st.button("Check selected variant in ClinVar and dbSNP",key="check_registry_variant"):
                try:
                    pos=codon_row["genomic_positions"][alternative["codon_base"]-1]
                    st.session_state["variant_registry_result"]={
                        "key":registry_variant_key,
                        "data":inspect_variant(loaded["annotation"]["seq_region_name"],pos,
                                               alternative["reference_base"],alternative["alternate_base"],
                                               loaded["annotation"]["strand"],lookup_assembly)
                    }
                except ValueError as exc:
                    st.error(str(exc))
            registry_result=st.session_state.get("variant_registry_result")
            if registry_result and registry_result.get("key")==registry_variant_key:
                data=registry_result["data"]
                st.markdown("#### ClinVar / dbSNP evidence for this exact genomic substitution")
                st.write("Genomic HGVS: "+data["variant"]["genomic_hgvs"])
                if data.get("curated_reference"):
                    known=data["curated_reference"]
                    st.success(f"Public reference record: {known['gene']} {known['coding_hgvs']} ({known['protein_short']}) | {known['rsid']}")
                    st.link_button("Review dbSNP reference",known["dbsnp_reference"])
                    st.link_button("Review ClinVar reference",known["clinvar_reference"])
                    st.caption(known["rflp_note"])
                st.write("dbSNP live-query status: "+data["dbsnp"]["status"])
                for rsid in data["dbsnp"].get("rsids",[]):
                    st.link_button(rsid,"https://www.ncbi.nlm.nih.gov/snp/"+rsid)
                st.write("ClinVar live-query status: "+data["clinvar"]["status"])
                for link in data["clinvar"].get("links",[]):
                    st.link_button("Review ClinVar candidate "+link.rstrip("/").split("/")[-1],link)
                st.caption(data["dbsnp"].get("note",""))
                st.caption(data["clinvar"].get("note",""))
                st.warning(data["novelty_assessment"])
                st.download_button("Download variant database evidence (JSON)",
                    __import__("json").dumps(data,indent=2),file_name="variant_registry_evidence.json",
                    mime="application/json")
            if st.button("Load reference genomic DNA for selected variant"):
                try:
                    with st.spinner("Retrieving reference flank and checking allele orientation..."):
                        template=fetch_refseq_assay(loaded["annotation"],loaded["cds"],codon_number,alternative["mutant_codon"],assembly=lookup_assembly) if lookup_provider.startswith("UCSC") else fetch_genomic_assay(loaded["annotation"],loaded["cds"],codon_number,alternative["mutant_codon"],assembly=lookup_assembly)
                    st.session_state["lookup_assay"]=template
                except ValueError as exc:
                    st.error(str(exc))
            template=st.session_state.get("lookup_assay")
            if template and template["metadata"]["transcript_id"]==tx["id"] and template["metadata"]["codon_number"]==codon_number and template["metadata"]["alternate_codon"]==alternative["mutant_codon"]:
                sequence=template["sequence"]
                derived_mutation=template["mutation"]
                variant_metadata={"variant":derived_mutation, "source":"Ensembl gene lookup",**template["metadata"]}
                st.success(f"Ready: {derived_mutation}, {len(sequence)} bp genomic flank, {template['metadata']['reference_codon']}→{template['metadata']['alternate_codon']}.")
                st.warning("This creates a reference-based candidate assay, not a genotype call or genome-wide primer specificity validation.")
elif mode=="WT marked motif → mutant motif":
    st.subheader("Mark the WT codon and enter its mutant replacement")
    marked_text=st.text_area("Wild-type DNA with exactly one marked region, e.g. AAA[GAG]CTG",height=170,
                            placeholder="Paste only DNA, with exactly one [WT motif], e.g. ACGT[GAG]TTAC...")
    alternate_motif=st.text_input("Mutant motif (same length; one base difference)","AAG")
    if marked_text.strip():
        try:
            sequence,derived_mutation,variant_metadata=parse_marked_reference(marked_text,alternate_motif)
            st.success(f"Calculated local substitution: {derived_mutation} (1-based, within supplied WT DNA).")
        except ValueError as exc:
            input_error=str(exc); st.error(input_error)
elif mode=="Complete WT and mutant sequences":
    st.subheader("Compare complete wild-type and mutant sequences")
    full_wt=st.text_area("Full wild-type sequence (A/C/G/T or one-record FASTA)",height=140)
    full_alt=st.text_area("Full mutant sequence (same length, exactly one changed base)",height=140)
    if full_wt.strip() and full_alt.strip():
        try:
            sequence,derived_mutation,variant_metadata=compare_wild_type_mutant(full_wt,full_alt)
            st.success(f"Derived variant: {derived_mutation}; both sequences have {len(sequence)} bp.")
        except ValueError as exc:
            input_error=str(exc); st.error(input_error)
elif mode=="Try example":
    st.subheader("Built-in example: TTN HMERF (200 bp)")
    st.write("Demonstration point mutation C101T in the included TTN exon sequence.")
    if SAMPLE.is_file():
        sequence=load_sequence(SAMPLE.read_text(encoding="utf-8"))
        st.code(sequence,language="text")
        st.info("Click Analyze example below to see an example report.")
    else:
        st.error("The example FASTA is missing from this checkout.")
elif mode=="Upload file":
    upload=st.file_uploader("Upload sequence or document",type=["fa","fasta","fna","txt","docx","pdf"])
    if upload:
        try:
            if upload.size > 10 * 1024 * 1024:
                raise ValueError("Upload exceeds the 10 MB limit.")
            candidates,alerts=extract_sequences(upload.name,upload.getvalue())
            for notice in alerts:
                st.info(notice)
            selected=st.selectbox("Select the correct reference sequence",
                                  list(range(len(candidates))),
                                  format_func=lambda i:candidates[i]["name"]+" ("+str(candidates[i]["length"])+" bp)")
            sequence=candidates[selected]["sequence"]
            st.success(f"Found {len(sequence):,} bp. Verify bases and reference coordinates.")
            st.text_area("Sequence preview (first 500 bp)",sequence[:500],height=110,disabled=True)
        except (ValueError,KeyError,RuntimeError) as exc:
            st.error("Cannot parse uploaded file: "+str(exc))
else:
    raw=st.text_area("Paste A/C/G/T sequence or single-record FASTA",height=160)
    if raw.strip():
        try:
            sequence=load_sequence(raw)
            st.caption(f"{len(sequence):,} bp of validated DNA")
        except ValueError as exc:
            st.error(str(exc))

# Never display a previous result for a different sequence or variant.
current_signature = (
    hashlib.sha256(sequence.encode("ascii")).hexdigest() if sequence else None,
    (derived_mutation if variant_metadata else mutation.strip().upper()) if input_error is None else "",
    str(variant_metadata) if variant_metadata else "",
    mode,
    variant_mode, reference_accession, genomic_start, genomic_strand, genome_assembly,
    name,
    max_enzymes,
)
if st.session_state.get("input_signature") != current_signature:
    st.session_state.pop("result", None)
    st.session_state.pop("sequence", None)
    st.session_state.pop("export", None)

if st.button("Analyze example" if mode=="Try example" else "Analyze sequence",
             type="primary",disabled=sequence is None or input_error is not None):
    try:
        assay_name=("TTN HMERF example" if name=="Variant" and mode=="Try example" else name)
        with st.spinner("Designing primers and screening enzyme patterns..."):
            metadata = None
            local_mutation = variant_metadata['variant'] if variant_metadata else mutation
            if variant_metadata is None and variant_mode.startswith("Genomic"):
                metadata = map_genomic_snv(sequence, mutation, reference_accession,
                                           int(genomic_start), genomic_strand, genome_assembly)
                local_mutation = metadata["local_mutation"]
            result=analyze_assay(sequence,local_mutation,assay_name,max_enzymes)
            if variant_metadata:
                result["variant_input"] = variant_metadata
                result["warnings"].append("Reference/variant input supplies a predicted substitution; no experimental genotype call is made.")
            if metadata:
                result["reference_mapping"] = metadata
                result["warnings"].append(metadata["mapping_warning"])
            result["enzymes"]=rank_enzymes(result["enzymes"])
        st.session_state["input_signature"]=current_signature
        st.session_state["result"]=result
        st.session_state["sequence"]=sequence
        st.session_state.pop("export",None)
        st.success("In-silico analysis completed. Review the ranking and validation limitations.")
    except (ValueError,RuntimeError,TypeError) as exc:
        st.session_state.pop("result",None)
        st.session_state.pop("export",None)
        st.error("Analysis could not complete: "+str(exc))

result=st.session_state.get("result")
if result and result.get("primers"):
    primer_options=result["primers"]
    primer_index=st.selectbox(
        "Choose primer pair for enzyme screening, virtual gel and reports",
        list(range(len(primer_options))),
        format_func=lambda i:f"Pair {i+1} — {primer_options[i]['product_size']} bp | F Tm {primer_options[i]['forward_tm']} °C | R Tm {primer_options[i]['reverse_tm']} °C",
        key="selected_primer_index"
    )
    if st.session_state.get("selected_primer_scope")!=(result["timestamp_utc"],primer_index):
        st.session_state["selected_primer_scope"]=(result["timestamp_utc"],primer_index)
        st.session_state["enzyme_index"]=0
        st.session_state.pop("export",None)
    candidate=primer_options[primer_index]
    original=st.session_state["sequence"]
    wt,mut=apply_mutation(original,result["position"]-1,result["ref"],result["alt"])
    left,right=candidate["forward_start"],candidate["reverse_start"]
    result=dict(result)
    result["selected_primer"]=candidate
    result["wt_amplicon"]=wt[left:right+1]
    result["mut_amplicon"]=mut[left:right+1]
    result["enzymes"]=rank_enzymes(find_differential_enzymes(
        result["wt_amplicon"],result["mut_amplicon"],max_results=max_enzymes))
    result["status"]="DESIGN_CANDIDATE" if result["enzymes"] else "NO_DIFFERENTIAL_ENZYME"
if result:
    st.divider()
    st.subheader(result["name"]+" — "+result["status"])
    m=st.columns(4)
    m[0].metric("Input DNA",f"{result['seq_length']:,} bp")
    m[1].metric("Variant",result["mutation"])
    m[2].metric("Primer pairs",len(result["primers"]))
    m[3].metric("Candidate enzymes",len(result["enzymes"]))
    meta=result.get("variant_input") or {}
    if meta.get("reference_codon"):
        st.info(f"Selected reference change: {meta['transcript_id']} | codon {meta['codon_number']} | {meta['reference_codon']} → {meta['alternate_codon']} | {meta['reference_aa']}{meta['codon_number']}{meta['alternate_aa']} | {meta['assembly']} {meta['genomic_chromosome']}:{meta['genomic_position']}")
    pages=st.tabs(["Primers & highlighted DNA","Restriction enzyme grades",
                   "Virtual gel & ladder","Download complete reports","Limitations"])
    with pages[0]:
        pair=result["selected_primer"]
        if pair:
            st.write("**Selected primer pair**")
            st.code("Forward 5'-"+pair["forward"]+"-3'\nReverse 5'-"+pair["reverse"]+"-3'")
            st.write(f"Predicted product: {pair['product_size']} bp | "
                     f"Forward Tm: {pair['forward_tm']} °C | Reverse Tm: {pair['reverse_tm']} °C")
            st.dataframe(pd.DataFrame(result["primers"]),hide_index=True,use_container_width=True)
            st.write("**Primer3 secondary-structure and local binding checks**")
            if pair.get("quality_metrics"):
                st.json(pair["quality_metrics"])
            st.caption("These checks are limited to the submitted DNA sequence, not the entire genome.")
            with tempfile.TemporaryDirectory() as tmp:
                output=render_sequence_figure(result,str(Path(tmp)/"primer_map.png"))
                st.image(Path(output).read_bytes(),caption="Primer-binding regions and local mutation coordinate")
            st.write("**Color-highlighted DNA reference and primer-binding sites**")
            components.html(sequence_markup(result),height=420,scrolling=True)
        else:
            st.error("No primer pair identified. No PCR-RFLP assay recommended.")
    with pages[1]:
        enzymes=result["enzymes"]
        if enzymes:
            frame=pd.DataFrame(enzymes)
            st.dataframe(frame[["priority_rank","enzyme","site","design_grade","priority_score",
                               "fragments_wt","fragments_mut","design_notes"]],
                         hide_index=True,use_container_width=True)
            top=enzymes[0]
            st.info(f"First computational choice: {top['enzyme']} — grade {top['design_grade']}, "
                    f"{top['priority_score']}/100. This is NOT a validated enzyme-activity ranking.")
            st.caption("A high computational score does NOT establish a best laboratory enzyme. Published assay validation and manufacturer documentation take precedence.")
            st.caption("Priority favors allele-discriminating, reasonably visible, simple cut patterns; "
                       "methylation sensitivity, star activity, reagent availability and price remain unverified.")
        else:
            st.info("No enzyme differentially cleaved the chosen primer amplicon.")
    with pages[2]:
        enzymes=result["enzymes"]
        options=list(range(len(enzymes)))
        enzyme_index=st.selectbox("Choose restriction enzyme",options,
            format_func=lambda i:(f"{enzymes[i]['priority_rank']}. {enzymes[i]['enzyme']} "
                                  f"— grade {enzymes[i]['design_grade']} ({enzymes[i]['priority_score']}/100)"),
            key="enzyme_index") if options else None
        current_enzyme=enzymes[enzyme_index] if enzyme_index is not None else None
        ladder_options=recommended_ladders(current_enzyme["fragments_wt"]+current_enzyme["fragments_mut"]) if current_enzyme else list(LADDERS)
        ladder_name=st.selectbox("Choose DNA molecular-weight ladder",ladder_options,key="ladder_name")
        ladder=LADDERS[ladder_name]
        st.caption(f"{ladder['vendor']} | {ladder['catalog']} | "
                   f"marker range {min(ladder['sizes'])}–{max(ladder['sizes'])} bp")
        st.caption("Marker sizes (bp): "+", ".join(map(str,ladder["sizes"])))
        if ladder["source"]:
            st.link_button("Ladder manufacturer's product information",ladder["source"])
        enzyme=enzymes[enzyme_index] if enzyme_index is not None else None
        if enzyme and result["selected_primer"]:
            with tempfile.TemporaryDirectory() as tmp:
                output=render_ladder_gel(result["selected_primer"]["product_size"],
                                         enzyme,ladder_name,str(Path(tmp)/"gel.png"))
                st.image(Path(output).read_bytes(),caption="Predicted complete-digestion gel")
            st.dataframe(pd.DataFrame([
                {"Genotype":"WT / WT","Predicted bands (bp)":", ".join(map(str,enzyme["fragments_wt"]))},
                {"Genotype":"WT / mutant","Predicted bands (bp)":", ".join(map(str,sorted(
                    enzyme["fragments_wt"]+enzyme["fragments_mut"],reverse=True)))},
                {"Genotype":"Mutant / mutant","Predicted bands (bp)":", ".join(map(str,enzyme["fragments_mut"]))}
            ]),hide_index=True,use_container_width=True)
            st.warning(ladder_notes(ladder_name,enzyme["fragments_wt"]+enzyme["fragments_mut"]))
        else:
            st.info("No enzyme-specific digest is available.")
    with pages[3]:
        st.markdown("**Complete downloadable package**")
        st.write("One ZIP containing PDF, Word, HTML, annotated sequence, primer map, "
                 "virtual gel, primers CSV, enzyme grades CSV, fragments CSV, JSON and FASTA.")
        ix=st.session_state.get("enzyme_index",0)
        selected_enzyme=result["enzymes"][ix] if result["enzymes"] and ix < len(result["enzymes"]) else None
        selected_ladder=st.session_state.get("ladder_name",next(iter(LADDERS)))
        key=(result["timestamp_utc"],st.session_state.get("selected_primer_index",0),selected_enzyme["enzyme"] if selected_enzyme else "",selected_ladder, str(result.get("variant_input",{})))
        if st.button("Prepare all reports",key="prepare_reports"):
            try:
                with st.spinner("Creating all downloadable results..."):
                    zipped,individual=build_export_bundle(result,st.session_state["sequence"],
                                                           selected_ladder,selected_enzyme)
                st.session_state["export"]=(key,zipped,individual)
                st.success("Reports prepared.")
            except Exception as exc:
                st.session_state.pop("export",None)
                st.error("Report generation failed: "+str(exc))
        export=st.session_state.get("export")
        if export and export[0]==key:
            st.download_button("Download ALL results (ZIP)",export[1],
                               file_name="genotype_complete_results.zip",
                               mime="application/zip",type="primary")
            for filename,body in export[2].items():
                st.download_button("Download "+filename,body,file_name=filename,
                                   key="download_"+filename)
        else:
            st.caption("Select enzyme and ladder on the preceding tab, then click Prepare all reports.")
    with pages[4]:
        for warning in result["warnings"]:
            st.warning(warning)
        st.caption("Not implemented: true HGVS transcript/genomic mapping, genome-wide primer specificity, "
                   "experimental conditions, validated dCAPS mismatch designs, clinical genotype calling.")
