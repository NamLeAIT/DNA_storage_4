
from __future__ import annotations
import os, tempfile, uuid
import pandas as pd
import streamlit as st
from ngs_fastq_v3 import FASTQSimulationConfig, export_simulated_fastq, simulate_fastq_from_library
from ngs_decode_v5 import load_ngs_prep_result_from_json
from utils_bits_v2 import ensure_dir

def render_sequencing_fastq():
    st.title("🧪 Sequencing / FASTQ")
    st.caption("Generate sequencing-style FASTQ reads with optional noise, or load external FASTQ files for downstream recovery.")
    mode=st.radio("Mode",["Simulate Sequencing","Load External FASTQ"], horizontal=True)
    prep=None
    if mode=="Simulate Sequencing":
        src=st.radio("Fragment library source",["From NGS Prep","Upload manifest JSON"], horizontal=True)
        if src=="From NGS Prep":
            prep=st.session_state.get("pipeline_ngs_prep_result")
            if prep: st.success(f"Loaded the latest fragment library with {len(prep.fragments):,} fragments.")
            else: st.warning("No fragment library is available from NGS Prep.")
        else:
            manifest=st.file_uploader("Upload manifest JSON", type=["json"])
            if manifest is not None:
                tmp=os.path.join(tempfile.gettempdir(),"seq_manifest_uploads",str(uuid.uuid4())); os.makedirs(tmp, exist_ok=True)
                path=os.path.join(tmp, manifest.name); open(path,"wb").write(manifest.getvalue())
                prep=load_ngs_prep_result_from_json(path)
                st.success(f"Loaded a fragment library with {len(prep.fragments):,} fragments.")
        c1,c2,c3,c4=st.columns(4)
        with c1: coverage=st.number_input("Coverage",1,500,20,1)
        with c2: paired_end=st.checkbox("Paired-end", value=True)
        with c3: read_length=st.number_input("Read length",20,300,150,1)
        with c4: mean_q=st.number_input("Mean Q",20,40,34,1)
        with st.expander("Show sequencing noise settings", expanded=False):
            e1,e2,e3,e4=st.columns(4)
            with e1: sub_rate=st.number_input("Substitution rate",0.0,0.05,0.003,0.001, format="%.3f")
            with e2: ins_rate=st.number_input("Insertion rate",0.0,0.02,0.0002,0.0001, format="%.4f")
            with e3: del_rate=st.number_input("Deletion rate",0.0,0.02,0.0002,0.0001, format="%.4f")
            with e4: seed=st.number_input("Random seed",0,999999,7,1)
        if st.button("Simulate FASTQ", type="primary", disabled=(prep is None), width="stretch"):
            sim_cfg=FASTQSimulationConfig(coverage=int(coverage), paired_end=bool(paired_end), read_length_nt=int(read_length), mean_q=int(mean_q), sub_rate=float(sub_rate), ins_rate=float(ins_rate), del_rate=float(del_rate), seed=int(seed))
            sim=simulate_fastq_from_library(prep, sim_cfg)
            out_dir=os.path.join("jobs","sequencing",str(uuid.uuid4())); ensure_dir(out_dir)
            sim_exports=export_simulated_fastq(sim, out_dir=out_dir, prefix="noisy_reads")
            st.session_state["pipeline_fastq_sim_result"]=sim; st.session_state["pipeline_fastq_exports"]=sim_exports
            st.success("Sequencing-style FASTQ files with simulated errors were generated successfully.")
    else:
        with st.container(border=True):
            r1=st.file_uploader("Upload R1 FASTQ(.gz)", type=["fastq","fq","gz"])
            r2=st.file_uploader("Upload R2 FASTQ(.gz) (optional)", type=["fastq","fq","gz"])
            st.caption("Use these files later in NGS Decode together with the manifest JSON exported by NGS Prep.")
            if st.button("Store FASTQ for NGS Decode", type="primary", disabled=(r1 is None), width="stretch"):
                out_dir=os.path.join("jobs","sequencing",str(uuid.uuid4())); ensure_dir(out_dir); exports={}
                if r1 is not None:
                    p=os.path.join(out_dir,r1.name); open(p,"wb").write(r1.getvalue()); exports["r1_fastq_gz"]=p
                if r2 is not None:
                    p=os.path.join(out_dir,r2.name); open(p,"wb").write(r2.getvalue()); exports["r2_fastq_gz"]=p
                st.session_state["pipeline_fastq_exports"]=exports; st.success("External FASTQ files were stored for downstream recovery.")
    sim=st.session_state.get("pipeline_fastq_sim_result"); exports=st.session_state.get("pipeline_fastq_exports")
    if sim:
        with st.container(border=True):
            st.markdown("### 📊 Injected Error Summary")
            s=sim.get("summary") or {}
            c1,c2,c3,c4,c5=st.columns(5)
            with c1: st.metric("Total Reads", f"{int(s.get('total_reads') or 0):,}")
            with c2: st.metric("Reads With Errors", f"{int(s.get('reads_with_errors') or 0):,}")
            with c3: st.metric("Substitutions", f"{int(s.get('total_substitutions') or 0):,}")
            with c4: st.metric("Insertions", f"{int(s.get('total_insertions') or 0):,}")
            with c5: st.metric("Deletions", f"{int(s.get('total_deletions') or 0):,}")
            ferr = s.get("fragment_error_summary") or []
            if ferr:
                df=pd.DataFrame(ferr)
                if "fragment_index" in df.columns:
                    df.insert(0,"Fragment", df["fragment_index"].map(lambda x:f"Fragment {int(x)}"))
                if "fragment_name" in df.columns: df=df.drop(columns=["fragment_name"])
                df=df.rename(columns={"fragment_index":"Number","is_parity_fragment":"Parity","reads_generated":"Reads Generated","reads_with_errors":"Reads With Errors","error_fraction":"Error Fraction","substitutions":"Substitutions","insertions":"Insertions","deletions":"Deletions"})
                st.dataframe(df, width="stretch", hide_index=True)
    if exports:
        with st.container(border=True):
            st.markdown("### 📤 FASTQ Output")
            if "r1_fastq_gz" in exports and os.path.exists(exports["r1_fastq_gz"]):
                st.download_button("Download R1 FASTQ.gz", data=open(exports["r1_fastq_gz"],"rb").read(), file_name=os.path.basename(exports["r1_fastq_gz"]), width="stretch")
            if "r2_fastq_gz" in exports and os.path.exists(exports["r2_fastq_gz"]):
                st.download_button("Download R2 FASTQ.gz", data=open(exports["r2_fastq_gz"],"rb").read(), file_name=os.path.basename(exports["r2_fastq_gz"]), width="stretch")
            if "summary_json" in exports and os.path.exists(exports["summary_json"]):
                st.download_button("Download FASTQ summary JSON", data=open(exports["summary_json"],"rb").read(), file_name=os.path.basename(exports["summary_json"]), width="stretch")
