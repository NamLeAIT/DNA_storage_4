
from __future__ import annotations
import os, uuid
from dataclasses import asdict
import pandas as pd, streamlit as st
from ngs_prep_v5 import NGSPrepConfig, PrimerConfig, clean_dna, gc_ratio, max_homopolymer, prepare_ngs_library, shannon_entropy, export_ngs_library, reconstruct_payload_from_fragments
from utils_bits_v2 import ensure_dir

def _load_latest_dna_from_history() -> str:
    history=st.session_state.get("history",[])
    if not history: return ""
    latest=history[-1]
    dna_path=((latest.get("artifacts",{}) or {}).get("dna_txt"))
    if dna_path and os.path.exists(dna_path):
        try: return clean_dna(open(dna_path,"r",encoding="utf-8",errors="ignore").read())
        except: return ""
    return ""

def _stats_block(seq:str):
    c1,c2,c3,c4=st.columns(4)
    with c1: st.metric("Total Length", f"{len(seq):,} nt" if seq else "—")
    with c2: st.metric("GC Ratio", f"{gc_ratio(seq):.4f}" if seq else "—")
    with c3: st.metric("Max Homopolymer", f"{max_homopolymer(seq)} nt" if seq else "—")
    with c4: st.metric("Shannon Entropy", f"{shannon_entropy(seq):.4f}" if seq else "—")

def _overhead(result):
    if not result or not result.fragments: return None
    first_non_parity=next((x for x in result.fragments if not x.is_parity_fragment), result.fragments[0])
    payload_len=first_non_parity.payload_length_nt
    primer_len=len(first_non_parity.forward_primer)+len(first_non_parity.reverse_primer)
    barcode_len=len(first_non_parity.sample_barcode); id_len=len(first_non_parity.fragment_id); prot_len=len(first_non_parity.checksum)
    total_library_nt=sum(x.fragment_length_nt for x in result.fragments); input_nt=int(result.input_sequence_length_nt or 0)
    return {"payload_len":payload_len,"primer_len":primer_len,"barcode_len":barcode_len,"id_len":id_len,"protection_len":prot_len,"fragment_len":first_non_parity.fragment_length_nt,"struct_overhead":first_non_parity.fragment_length_nt-payload_len,"total_library_nt":total_library_nt,"input_nt":input_nt,"expansion":(total_library_nt/input_nt if input_nt else None)}

def render_ngs_prep():
    st.title("🧫 NGS Prep")
    st.caption("Turn a DNA payload into a sequencing-ready fragment library.")
    with st.container(border=True):
        st.markdown("### What happens in this step?")
        st.markdown("- The DNA payload is split into sequencing fragments.\n- Each fragment receives primers, a fragment identifier, and the selected protection method.\n- The app exports a manifest and fragment library for sequencing and recovery.")
    left,right=st.columns([1,1], gap="small")
    with left:
        with st.container(border=True):
            st.markdown("### 📥 Input")
            source=st.radio("DNA Source",["From Encode tab","Upload DNA payload file","Paste DNA directly"], horizontal=True)
            seq=""
            if source=="From Encode tab":
                seq=_load_latest_dna_from_history()
                if seq: st.success(f"Loaded DNA payload from Encode ({len(seq):,} nt).")
                else: st.warning("No DNA payload is available from the current Encode history.")
            elif source=="Upload DNA payload file":
                up=st.file_uploader("Upload DNA payload file", type=["txt","dna","fasta","fa"])
                if up is not None:
                    seq=clean_dna(up.getvalue().decode("utf-8",errors="ignore"))
                    if seq: st.success(f"Loaded DNA payload ({len(seq):,} nt).")
            else:
                direct=st.text_area("Paste DNA sequence", height=180)
                seq=clean_dna(direct)
                if seq: st.success(f"Loaded DNA payload ({len(seq):,} nt).")
        with st.container(border=True):
            st.markdown("### 📊 Input DNA Statistics"); _stats_block(seq)
    with right:
        with st.container(border=True):
            st.markdown("### ⚙️ Core Settings")
            c1,c2=st.columns(2)
            with c1:
                payload_len=st.number_input("Payload length per fragment",20,200,80,1)
                protection_mode=st.selectbox("Protection mode", ["None","CRC","Parity Fragments","Reed-Solomon"], index=0)
            with c2:
                primer_source=st.selectbox("Primer source", ["Default (Illumina-compatible)","Custom"], index=0)
                last_policy=st.selectbox("Last fragment handling", ["Keep shorter","Pad to full length"], index=0)
            with st.expander("Show advanced settings", expanded=False):
                a1,a2=st.columns(2)
                with a1:
                    fragment_id_len=st.number_input("Fragment ID length",0,24,8,1)
                    checksum_len=st.number_input("Protection length (nt)",0,16,8 if protection_mode=="CRC" else 0,1)
                    parity_group_size=st.number_input("Parity group size",2,16,4,1, disabled=(protection_mode!="Parity Fragments"))
                    rs_data_shards=st.number_input("RS data shards",2,12,4,1, disabled=(protection_mode!="Reed-Solomon"))
                    rs_parity_shards=st.number_input("RS parity shards",1,6,2,1, disabled=(protection_mode!="Reed-Solomon"))
                with a2:
                    sample_barcode=st.text_input("Sample barcode", value="")
                    target_gc=st.slider("Target GC",0.35,0.65,0.50,0.01)
                if primer_source=="Default (Illumina-compatible)":
                    fwd="ACACGACGCTCTTCCGATCT"; rev="AGATCGGAAGAGCACACGTCT"
                else:
                    fwd=st.text_input("Forward primer (5')", value="ACACGACGCTCTTCCGATCT")
                    rev=st.text_input("Reverse primer (3')", value="AGATCGGAAGAGCACACGTCT")
            if primer_source=="Default (Illumina-compatible)":
                fwd="ACACGACGCTCTTCCGATCT"; rev="AGATCGGAAGAGCACACGTCT"
            st.caption("Output of this step is a sequencing-ready fragment library. Use the manifest JSON and fragment FASTA in the next steps.")
    if st.button("Generate NGS Fragments", type="primary", disabled=(not seq), width="stretch"):
        cfg=NGSPrepConfig(
            payload_length_nt=int(payload_len),
            fragment_id_length_nt=int(fragment_id_len),
            sample_barcode=sample_barcode,
            checksum_length_nt=(int(checksum_len) if protection_mode=="CRC" else 0),
            primer=PrimerConfig(forward=fwd, reverse=rev),
            last_fragment_policy=("pad_to_full" if last_policy=="Pad to full length" else "keep_shorter"),
            target_gc=float(target_gc),
            protection_mode=protection_mode,
            parity_group_size=int(parity_group_size),
            rs_data_shards=int(rs_data_shards),
            rs_parity_shards=int(rs_parity_shards),
        )
        result=prepare_ngs_library(seq, cfg)
        out_dir=os.path.join("jobs","ngs",str(uuid.uuid4())); ensure_dir(out_dir)
        exports=export_ngs_library(result,out_dir,prefix="ngs_library")
        st.session_state["pipeline_ngs_prep_result"]=result
        st.session_state["pipeline_ngs_prep_exports"]=exports
        st.session_state["pipeline_payload_before_fastq"]=reconstruct_payload_from_fragments(result.fragments)
        st.session_state["pipeline_protection_mode"]=protection_mode
        st.success(f"Generated {result.summary['total_fragments']:,} sequencing fragments.")
    result=st.session_state.get("pipeline_ngs_prep_result"); exports=st.session_state.get("pipeline_ngs_prep_exports")
    if result:
        with st.container(border=True):
            st.markdown("### 📤 Fragment Library Summary")
            s=result.summary; c1,c2,c3,c4=st.columns(4)
            with c1: st.metric("Total Fragments", f"{s['total_fragments']:,}")
            with c2: st.metric("Data Fragments", f"{int(s.get('data_fragments') or 0):,}")
            with c3: st.metric("Protection Fragments", f"{int(s.get('parity_fragments') or 0):,}")
            with c4: st.metric("Protection mode", st.session_state.get("pipeline_protection_mode","None"))
            ov=_overhead(result)
            if ov:
                st.markdown("#### Protection Overhead and Storage Efficiency")
                o1,o2,o3,o4,o5=st.columns(5)
                with o1: st.metric("Payload per fragment", f"{ov['payload_len']} nt")
                with o2: st.metric("Structural overhead", f"{ov['struct_overhead']} nt")
                with o3: st.metric("Protection length", f"{ov['protection_len']} nt")
                with o4: st.metric("Total library length", f"{ov['total_library_nt']:,} nt")
                with o5: st.metric("Expansion factor", f"{ov['expansion']:.2f}×" if ov['expansion'] else "—")
            df=pd.DataFrame([asdict(x) for x in result.fragments])
            view=df[["fragment_index","is_parity_fragment","payload_length_nt","fragment_length_nt","gc_ratio","max_homopolymer","score","difficulty"]].copy()
            view.insert(0,"Fragment", view["fragment_index"].map(lambda x:f"Fragment {int(x)}"))
            view=view.rename(columns={"fragment_index":"Number","is_parity_fragment":"Parity","payload_length_nt":"Payload Length (nt)","fragment_length_nt":"Fragment Length (nt)","gc_ratio":"GC Ratio","max_homopolymer":"Max Homopolymer","score":"Score","difficulty":"Difficulty"})
            st.dataframe(view, width="stretch", hide_index=True)
