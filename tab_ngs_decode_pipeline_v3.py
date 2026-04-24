
from __future__ import annotations
import os, tempfile, uuid
import pandas as pd, streamlit as st
from ngs_decode_v6 import decode_dna_payload_to_file, dna_stats, load_ngs_prep_result_from_json, precheck_decode_payload, reconstruct_payload_exact, recover_payload_from_fastq
from ngs_prep_v5 import clean_dna
from preview_utils_v1 import render_file_preview

def _stats_block(seq:str):
    s=dna_stats(seq); c1,c2,c3,c4=st.columns(4)
    with c1: st.metric("Total Length", f"{s['length_nt']:,} nt" if seq else "—")
    with c2: st.metric("GC Ratio", f"{s['gc_ratio']:.4f}" if seq else "—")
    with c3: st.metric("Max Homopolymer", f"{s['max_homopolymer']} nt" if seq else "—")
    with c4: st.metric("Shannon Entropy", f"{s['shannon_entropy']:.4f}" if seq else "—")

def render_ngs_decode():
    st.title("📥 NGS Decode")
    st.caption("Reconstruct DNA from a clean payload, fragment library, or FASTQ reads, then optionally decode the recovered DNA back into a file.")
    left,right=st.columns([1,1], gap="small")
    with left:
        with st.container(border=True):
            st.markdown("### 📥 Decode Source")
            source=st.radio("Source",["DNA payload before sequencing","Fragment library before FASTQ","FASTQ after sequencing or simulation","Recovered DNA text file"])
            prep=None; fastq_r1_path=None; fastq_r2_path=None; payload=""
            if source=="DNA payload before sequencing":
                payload=clean_dna(st.session_state.get("pipeline_payload_before_fastq") or "")
                if payload: st.success(f"Loaded the pre-sequencing DNA payload ({len(payload):,} nt).")
                else:
                    up=st.file_uploader("Upload DNA payload text file", type=["txt","fasta","fa","dna"])
                    if up is not None:
                        payload=clean_dna(up.getvalue().decode("utf-8",errors="ignore"))
                        if payload: st.success(f"Loaded the DNA payload ({len(payload):,} nt).")
            elif source=="Fragment library before FASTQ":
                st.caption("Required file: manifest JSON from NGS Prep.")
                use_latest=st.checkbox("Use the latest NGS Prep result from this session", value=True)
                if use_latest and st.session_state.get("pipeline_ngs_prep_result") is not None:
                    prep=st.session_state.get("pipeline_ngs_prep_result"); st.success(f"Loaded the latest fragment library with {len(prep.fragments):,} fragments.")
                else:
                    manifest=st.file_uploader("Upload manifest JSON", type=["json"])
                    if manifest is not None:
                        tmp=os.path.join(tempfile.gettempdir(),"dec_manifest_only",str(uuid.uuid4())); os.makedirs(tmp,exist_ok=True)
                        path=os.path.join(tmp,manifest.name); open(path,"wb").write(manifest.getvalue())
                        prep=load_ngs_prep_result_from_json(path); st.success(f"Loaded the fragment library manifest with {len(prep.fragments):,} fragments.")
                if prep is not None: payload=reconstruct_payload_exact(prep)
            elif source=="FASTQ after sequencing or simulation":
                st.caption("Required files: manifest JSON from NGS Prep, one R1 FASTQ file, and optionally one R2 FASTQ file.")
                use_latest=st.checkbox("Use the latest sequencing run from this session", value=True)
                if use_latest and st.session_state.get("pipeline_ngs_prep_result") is not None:
                    prep=st.session_state.get("pipeline_ngs_prep_result")
                    ex=st.session_state.get("pipeline_fastq_exports") or {}
                    fastq_r1_path=ex.get("r1_fastq_gz"); fastq_r2_path=ex.get("r2_fastq_gz")
                    if prep: st.success(f"Loaded the latest fragment library with {len(prep.fragments):,} fragments.")
                    if fastq_r1_path and os.path.exists(fastq_r1_path): st.info(f"Loaded R1 FASTQ: {os.path.basename(fastq_r1_path)}")
                    else: st.warning("No R1 FASTQ file is available from the current session.")
                else:
                    manifest=st.file_uploader("Upload manifest JSON", type=["json"], key="dec_fastq_manifest_v2")
                    r1=st.file_uploader("Upload R1 FASTQ(.gz)", type=["fastq","fq","gz"], key="dec_fastq_r1_v2")
                    r2=st.file_uploader("Upload R2 FASTQ(.gz) (optional)", type=["fastq","fq","gz"], key="dec_fastq_r2_v2")
                    if manifest is not None:
                        tmp=os.path.join(tempfile.gettempdir(),"dec_fastq_files",str(uuid.uuid4())); os.makedirs(tmp,exist_ok=True)
                        mpath=os.path.join(tmp,manifest.name); open(mpath,"wb").write(manifest.getvalue())
                        prep=load_ngs_prep_result_from_json(mpath)
                        if r1 is not None:
                            fastq_r1_path=os.path.join(tmp,r1.name); open(fastq_r1_path,"wb").write(r1.getvalue())
                        if r2 is not None:
                            fastq_r2_path=os.path.join(tmp,r2.name); open(fastq_r2_path,"wb").write(r2.getvalue())
                        st.success(f"Loaded the fragment library manifest with {len(prep.fragments):,} fragments.")
            else:
                up=st.file_uploader("Upload recovered DNA text file", type=["txt","fasta","fa","dna"])
                if up is not None:
                    payload=clean_dna(up.getvalue().decode("utf-8",errors="ignore"))
                    if payload: st.success(f"Loaded recovered DNA text ({len(payload):,} nt).")
        with st.container(border=True):
            st.markdown("### ⚙️ Recovery Method")
            available_methods=["Consensus reconstruction"]
            prot=st.session_state.get("pipeline_protection_mode") or ((prep.config.get("protection_mode") if prep else None) or "None")
            if prot in ("CRC","Parity Fragments","Reed-Solomon"):
                available_methods += ["Consensus + checksum validation","ECC-assisted recovery"]
            method=st.selectbox("Method", available_methods, index=0)
            with st.expander("Show advanced decode settings", expanded=False):
                mode_label=st.selectbox("DNA mapping", ["Simple Mapping","Rule-based"], index=0)
                mode_codec="SIMPLE" if mode_label=="Simple Mapping" else "TABLE"
                scheme_name="RINF_B16"; init_dimer="TA"; whiten=False
                if mode_codec=="TABLE":
                    scheme_name=st.selectbox("Rule", ["R0_B9","R1_B12","R2_B15","RINF_B16"], index=3)
                    init_dimer=st.selectbox("Initial dimer", ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"], index=9)
                preferred_stem=st.text_input("Recovered output name", value="recovered_from_ngs")
            # defaults hidden
            if "mode_codec" not in locals():
                mode_codec="SIMPLE"; scheme_name="RINF_B16"; init_dimer="TA"; whiten=False; preferred_stem="recovered_from_ngs"
            r1,r2=st.columns(2)
            with r1: run_recovery=st.button("Run Recovery", type="primary", width="stretch")
            with r2: run_decode=st.button("Decode Recovered Payload", width="stretch")
    with right:
        with st.container(border=True):
            st.markdown("### 📊 Selected Input"); _stats_block(payload)
    if run_recovery:
        report=None
        if source in ("DNA payload before sequencing","Recovered DNA text file"):
            if not payload:
                st.error("No DNA payload is available for recovery."); return
            report={"mode":"payload_direct","total_reads":0,"total_fragments":0,"recovered_fragments":0,"recovery_fraction":1.0,"rows":[],"protection_mode":"None"}
        elif source=="Fragment library before FASTQ":
            if prep is None:
                st.error("Please provide a valid fragment library manifest."); return
            payload=reconstruct_payload_exact(prep)
            report={"mode":"fragment_library_direct","total_reads":0,"total_fragments":len([f for f in prep.fragments if not f.is_parity_fragment]),"recovered_fragments":len([f for f in prep.fragments if not f.is_parity_fragment]),"recovery_fraction":1.0,"rows":[],"protection_mode":(prep.config or {}).get("protection_mode","None")}
        else:
            if prep is None or not fastq_r1_path:
                st.error("Please provide the NGS manifest and at least one R1 FASTQ file."); return
            payload, report = recover_payload_from_fastq(prep, fastq_r1_path, fastq_r2_path, method=method)
        st.session_state["pipeline_decode_payload"]=payload
        st.session_state["pipeline_decode_report"]=report
        st.session_state["pipeline_decode_source"]=source
        st.session_state["pipeline_decode_method"]=method
        st.session_state["pipeline_decode_protection_mode"]=report.get("protection_mode") if isinstance(report, dict) else "None"
        pre=precheck_decode_payload(payload, scheme_name=scheme_name, mode_codec=mode_codec, seed="rn", init_dimer=init_dimer, whiten=whiten, remove_leading_one=True, target_gc=0.50, w_gc=0.0, w_motif=0.0, ks=(4,6))
        st.session_state["pipeline_decode_precheck"]=pre
    if run_decode:
        payload=clean_dna(st.session_state.get("pipeline_decode_payload") or "")
        if not payload:
            st.error("Run Recovery first to generate a payload."); return
        stats, restored=decode_dna_payload_to_file(payload, scheme_name=scheme_name, mode_codec=mode_codec, seed="rn", init_dimer=init_dimer, whiten=whiten, remove_leading_one=True, target_gc=0.50, w_gc=0.0, w_motif=0.0, ks=(4,6), preferred_stem=preferred_stem)
        st.session_state["pipeline_decode_final"]=stats; st.session_state["pipeline_decode_restored"]=restored
    report=st.session_state.get("pipeline_decode_report"); payload=clean_dna(st.session_state.get("pipeline_decode_payload") or ""); precheck=st.session_state.get("pipeline_decode_precheck"); final_stats=st.session_state.get("pipeline_decode_final"); restored=st.session_state.get("pipeline_decode_restored")
    if report:
        with st.container(border=True):
            st.markdown("### 📈 Recovery Summary")
            c1,c2,c3,c4,c5=st.columns(5)
            with c1: st.metric("Protection mode", str(report.get("protection_mode") or "None"))
            with c2: st.metric("Total Reads", f"{int(report.get('total_reads') or 0):,}")
            with c3: st.metric("Expected Fragments", f"{int(report.get('total_fragments') or 0):,}")
            with c4: st.metric("Recovered Fragments", f"{int(report.get('recovered_fragments') or 0):,}")
            with c5: st.metric("Recovery Rate", f"{float(report.get('recovery_fraction') or 0.0)*100:.1f}%")
            rows=report.get("rows") or []
            if rows:
                df=pd.DataFrame(rows)
                if "fragment_index" in df.columns: df.insert(0,"Fragment", df["fragment_index"].map(lambda x:f"Fragment {int(x)}"))
                if "fragment_name" in df.columns: df=df.drop(columns=["fragment_name"])
                df=df.rename(columns={"fragment_index":"Number","is_parity_fragment":"Parity","supporting_positions":"Supporting Positions","difficulty":"Difficulty","score":"Score","payload_length_nt":"Payload Length (nt)","recovered":"Recovered","status":"Status","mismatch_fraction_vs_manifest":"Mismatch Fraction","crc_ok":"CRC OK"})
                st.dataframe(df, width="stretch", hide_index=True)
    if payload:
        with st.container(border=True):
            st.markdown("### 🧬 Recovered DNA Payload"); _stats_block(payload)
            st.text_area("Recovered DNA preview", value=payload[:12000], height=240)
            st.download_button("Download recovered DNA payload", data=(payload+"\n").encode("utf-8"), file_name="recovered_dna_payload.txt", width="stretch")
    if precheck:
        with st.container(border=True):
            st.markdown("### 🧪 Decode Check")
            if precheck.get("ok"): st.success("The reconstructed DNA passes the decode and zlib integrity check.")
            else: st.error("The reconstructed DNA does not pass the decode and zlib integrity check.")
            st.json(precheck, expanded=False)
    if final_stats:
        with st.container(border=True):
            st.markdown("### 📥 Final Decode Result")
            d1,d2,d3=st.columns(3)
            with d1: st.metric("Decoded Bits", f"{int(final_stats.get('decoded_bits_len') or 0):,}")
            with d2: st.metric("Inner Bytes", f"{int(final_stats.get('inner_bytes_len') or 0):,}")
            with d3:
                z=(final_stats.get("zlib") or {})
                st.metric("zlib Integrity", "OK" if z.get("integrity_ok") else "Failed")
            st.json(final_stats, expanded=False)
            if restored and os.path.exists(restored):
                st.success(f"Recovered file: {os.path.basename(restored)}")
                with open(restored,"rb") as f: st.download_button("Download recovered file", data=f.read(), file_name=os.path.basename(restored), width="stretch")
                render_file_preview(restored, "Recovered File Preview")
