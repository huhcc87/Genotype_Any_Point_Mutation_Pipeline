"""Guided PCR-RFLP design UI: document upload, ranked enzymes, gel, complete reports."""
from pathlib import Path
import tempfile
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from src.assay import analyze_assay
from src.document_input import extract_sequences
from src.sequence import load_sequence
from src.ranking import rank_enzymes
from src.ladders import LADDERS,ladder_notes
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
    mode=st.radio("Input source",["Try example","Upload file","Paste sequence"])
    name=st.text_input("Assay name","Variant")
    mutation=st.text_input("Mutation (local sequence position)","C101T")
    max_enzymes=st.slider("Maximum restriction enzymes",3,50,20)
    st.divider()
    st.caption("Accepted: .fa, .fasta, .fna, .txt, .docx, .pdf (10 MB maximum).")
    st.caption("Do not upload identifiable patient data to any publicly hosted application.")

sequence=None
if mode=="Try example":
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

if st.button("Analyze example" if mode=="Try example" else "Analyze sequence",
             type="primary",disabled=sequence is None):
    try:
        assay_name=("TTN HMERF example" if name=="Variant" and mode=="Try example" else name)
        with st.spinner("Designing primers and screening enzyme patterns..."):
            result=analyze_assay(sequence,mutation,assay_name,max_enzymes)
            result["enzymes"]=rank_enzymes(result["enzymes"])
        st.session_state["result"]=result
        st.session_state["sequence"]=sequence
        st.session_state.pop("export",None)
        st.success("In-silico analysis completed. Review the ranking and validation limitations.")
    except (ValueError,RuntimeError,TypeError) as exc:
        st.session_state.pop("result",None)
        st.session_state.pop("export",None)
        st.error("Analysis could not complete: "+str(exc))

result=st.session_state.get("result")
if result:
    st.divider()
    st.subheader(result["name"]+" — "+result["status"])
    m=st.columns(4)
    m[0].metric("Input DNA",f"{result['seq_length']:,} bp")
    m[1].metric("Variant",result["mutation"])
    m[2].metric("Primer pairs",len(result["primers"]))
    m[3].metric("Candidate enzymes",len(result["enzymes"]))
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
        ladder_name=st.selectbox("Choose DNA molecular-weight ladder",list(LADDERS),key="ladder_name")
        ladder=LADDERS[ladder_name]
        st.caption(f"{ladder['vendor']} | {ladder['catalog']} | "
                   f"marker range {min(ladder['sizes'])}–{max(ladder['sizes'])} bp")
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
        key=(result["timestamp_utc"],selected_enzyme["enzyme"] if selected_enzyme else "",selected_ladder)
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
