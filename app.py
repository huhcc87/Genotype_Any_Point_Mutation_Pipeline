import io, json, tempfile, zipfile
from pathlib import Path
import pandas as pd
import streamlit as st
from src.assay import analyze_assay
from src.gel import draw_virtual_gel

st.set_page_config(page_title="Genotype PCR-RFLP", page_icon="🧬", layout="wide")
st.title("Genotype | PCR-RFLP assay designer")
st.warning("Research-use predictions only. Validate assay specificity, coordinate mapping, restriction conditions and expected bands experimentally.")
with st.sidebar:
    name=st.text_input("Assay name","Variant")
    variant=st.text_input("Mutation (relative to supplied sequence)","C101T")
    limit=st.slider("Enzymes per primer pair",1,30,10)
    run=st.button("Analyze",type="primary")
file=st.file_uploader("Optional single-record FASTA upload",type=["fa","fasta","fna","txt"])
sequence=st.text_area("DNA sequence or FASTA text",height=160,disabled=file is not None)
if run:
    try:
        contents=file.getvalue().decode("utf-8") if file else sequence
        st.session_state["result"]=analyze_assay(contents,variant,name,limit)
    except (ValueError,UnicodeError) as exc:
        st.session_state.pop("result",None)
        st.error(str(exc))
result=st.session_state.get("result")
if result:
    st.subheader(result["name"]+" — "+result["status"])
    a,b,c=st.columns(3)
    a.metric("Sequence (bp)",result["seq_length"])
    b.metric("Primer pairs",len(result["primers"]))
    c.metric("Differential enzymes",len(result["enzymes"]))
    t1,t2,t3,t4=st.tabs(["Primers","Restriction digests","Virtual gel","Downloads & limitations"])
    with t1:
        st.dataframe(pd.DataFrame(result["primers"]),use_container_width=True)
    with t2:
        st.dataframe(pd.DataFrame(result["enzymes"]),use_container_width=True)
    with t3:
        if result["enzymes"] and result["selected_primer"]:
            e=result["enzymes"][0]
            with tempfile.TemporaryDirectory() as tmp:
                path=draw_virtual_gel(result["selected_primer"]["product_size"],e["cuts_wt"],e["cuts_mut"],e["enzyme"],str(Path(tmp)/"gel.png"))
                img=Path(path).read_bytes()
            st.image(img,caption="Complete-digestion prediction; not observed laboratory bands")
            st.download_button("Download PNG",img,file_name="virtual_gel.png",mime="image/png")
        else:
            st.info("No discriminatory restriction digest for this amplicon.")
    with t4:
        for warning in result["warnings"]:
            st.warning(warning)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/"assay.json").write_text(json.dumps(result,indent=2))
            pd.DataFrame(result["enzymes"]).to_csv(root/"enzymes.csv",index=False)
            pd.DataFrame(result["primers"]).to_csv(root/"primers.csv",index=False)
            stream=io.BytesIO()
            with zipfile.ZipFile(stream,"w",zipfile.ZIP_DEFLATED) as archive:
                for file_path in root.iterdir():
                    archive.write(file_path,arcname=file_path.name)
        st.download_button("Download complete analysis ZIP",stream.getvalue(),file_name="genotype_analysis.zip",mime="application/zip")
