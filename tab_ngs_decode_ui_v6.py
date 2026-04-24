
from __future__ import annotations

import os
import tempfile
import uuid

import pandas as pd
import streamlit as st

from ngs_decode_v4 import (
    decode_dna_payload_to_file,
    dna_stats,
    load_ngs_prep_result_from_json,
    precheck_decode_payload,
    reconstruct_payload_exact,
    reconstruct_payload_from_fastq_standard,
    rescue_payload_from_fastq_advanced,
)
from ngs_prep_v3 import clean_dna


def _latest_bundle() -> dict:
    return st.session_state.get("ngs_latest_exports") or {}


def _latest_prep():
    return _latest_bundle().get("prep_result")


def _latest_fastq_exports():
    return _latest_bundle().get("sim_exports") or st.session_state.get("ngs_fastq_exports") or {}


def _latest_payload_before_fastq() -> str:
    return clean_dna((_latest_bundle().get("payload_before_fastq") or ""))


def _stats_block(seq: str):
    stats = dna_stats(seq)
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Total Length", f"{stats['length_nt']:,} nt" if seq else "—")
    with c2:
        st.metric("GC Ratio", f"{stats['gc_ratio']:.4f}" if seq else "—")
    with c3:
        st.metric("Max Homopolymer", f"{stats['max_homopolymer']} nt" if seq else "—")
    with c4:
        st.metric("Shannon Entropy", f"{stats['shannon_entropy']:.4f}" if seq else "—")


def render_ngs_decode():
    st.title("📥 NGS Decode")
    st.caption("Recover DNA from sequencing-style inputs, inspect the reconstructed payload, and optionally decode the payload back into a file.")

    left, right = st.columns([1.0, 1.0], gap="small")

    with left:
        with st.container(border=True):
            st.markdown("### 🧭 Recovery Level")
            level = st.radio(
                "Recovery method",
                ["Consensus recovery", "Manifest-guided rescue"],
                horizontal=True,
                key="ngs_dec_level_v4",
            )
            st.info(
                "Standard Test is for ordinary users: inspect fragment DNA, simulate sequencing noise, and see whether normal recovery works. "
                "Advanced Rescue is a stronger path that tries to rescue sequencing-affected DNA before final decode."
            )

        with st.container(border=True):
            st.markdown("### 📥 Decode Source")
            source = st.radio(
                "Source",
                [
                    "From DNA before FASTQ",
                    "From latest FASTQ run",
                    "Import FASTQ files",
                    "Import DNA text file",
                ],
                key="ngs_dec_source_v4",
            )

            prep = None
            fastq_r1_path = None
            fastq_r2_path = None
            selected_payload = ""

            if source == "From DNA before FASTQ":
                prep = _latest_prep()
                selected_payload = _latest_payload_before_fastq()
                if selected_payload:
                    st.success(f"Loaded DNA payload before FASTQ ({len(selected_payload):,} nt)")
                else:
                    st.warning("No DNA-before-FASTQ payload is available in the current session.")

            elif source == "From latest FASTQ run":
                prep = _latest_prep()
                sim_exports = _latest_fastq_exports()
                fastq_r1_path = sim_exports.get("r1_fastq_gz")
                fastq_r2_path = sim_exports.get("r2_fastq_gz")
                if prep:
                    st.success(f"Loaded the latest sequencing library with {len(prep.fragments):,} fragments.")
                else:
                    st.warning("No NGS Prep result found in the current session.")
                if fastq_r1_path and os.path.exists(fastq_r1_path):
                    st.info(f"Loaded the latest simulated R1 FASTQ file: {os.path.basename(fastq_r1_path)}")
                else:
                    st.warning("No FASTQ from latest NGS run is available.")

            elif source == "Import FASTQ files":
                manifest = st.file_uploader("Upload NGS manifest JSON", type=["json"], key="ngs_dec_manifest_v4")
                r1 = st.file_uploader("Upload R1 FASTQ(.gz)", type=["fastq", "fq", "gz"], key="ngs_dec_r1_v4")
                r2 = st.file_uploader("Upload R2 FASTQ(.gz) (optional)", type=["fastq", "fq", "gz"], key="ngs_dec_r2_v4")
                st.caption("Files required for this option: the manifest JSON exported by NGS Preparation, one R1 FASTQ file, and optionally one R2 FASTQ file.")

                if manifest is not None:
                    tmp_dir = os.path.join(tempfile.gettempdir(), "ngs_decode_v4_uploads", str(uuid.uuid4()))
                    os.makedirs(tmp_dir, exist_ok=True)
                    manifest_path = os.path.join(tmp_dir, manifest.name)
                    with open(manifest_path, "wb") as f:
                        f.write(manifest.getvalue())
                    prep = load_ngs_prep_result_from_json(manifest_path)
                    st.success(f"Loaded the library manifest with {len(prep.fragments):,} fragments.")
                    if r1 is not None:
                        fastq_r1_path = os.path.join(tmp_dir, r1.name)
                        with open(fastq_r1_path, "wb") as f:
                            f.write(r1.getvalue())
                    if r2 is not None:
                        fastq_r2_path = os.path.join(tmp_dir, r2.name)
                        with open(fastq_r2_path, "wb") as f:
                            f.write(r2.getvalue())

            else:
                up = st.file_uploader("Upload DNA text file", type=["txt", "fasta", "fa", "dna"], key="ngs_dec_dna_file_v4")
                if up is not None:
                    selected_payload = clean_dna(up.getvalue().decode("utf-8", errors="ignore"))
                    if selected_payload:
                        st.success(f"Loaded DNA file ({len(selected_payload):,} nt)")

        with st.container(border=True):
            st.markdown("### ⚙️ Decode Settings")
            mode_label = st.selectbox("Mode", ["Simple Mapping", "Rule-base"], index=0, key="ngs_dec_mode_v4")
            mode_codec = "SIMPLE" if mode_label == "Simple Mapping" else "TABLE"
            scheme_name = "RINF_B16"
            init_dimer = "TA"
            if mode_codec == "TABLE":
                scheme_name = st.selectbox("DNA rule", ["R0_B9", "R1_B12", "R2_B15", "RINF_B16"], index=3, key="ngs_dec_scheme_v4")
                init_dimer = st.selectbox("Initial dimer", ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"], index=9, key="ngs_dec_init_v4")
            preferred_stem = st.text_input("Recovered output name", value="recovered_from_ngs", key="ngs_dec_stem_v4")

            if level == "Consensus recovery":
                run = st.button("Recover DNA from reads", type="primary", width="stretch", key="run_std_v4")
            else:
                run = st.button("Recover DNA with rescue", type="primary", width="stretch", key="run_adv_v4")

    with right:
        with st.container(border=True):
            st.markdown("### 📊 Selected Input")
            preview = ""
            if source == "From DNA before FASTQ":
                preview = selected_payload
            elif source in ("From latest FASTQ run", "Import FASTQ files") and prep is not None:
                preview = reconstruct_payload_exact(prep)
            elif source == "Import DNA text file":
                preview = selected_payload
            _stats_block(preview)

    if run:
        report = None
        payload = ""
        precheck = None
        final_stats = None
        restored = None

        if source == "From DNA before FASTQ":
            payload = selected_payload
            if not payload:
                st.error("No pre-sequencing DNA payload is available in the current session.")
                return

        elif source in ("From latest FASTQ run", "Import FASTQ files"):
            if prep is None or not fastq_r1_path:
                st.error("Please provide the NGS library manifest and at least one R1 FASTQ file.")
                return
            if level == "Consensus recovery":
                payload, report = reconstruct_payload_from_fastq_standard(prep, fastq_r1_path, fastq_r2_path)
            else:
                payload, report = rescue_payload_from_fastq_advanced(prep, fastq_r1_path, fastq_r2_path)

        else:
            payload = selected_payload
            if not payload:
                st.error("Please upload a DNA text file.")
                return

        st.session_state["ngs_dec_payload_v4"] = payload
        st.session_state["ngs_dec_report_v4"] = report
        st.session_state["ngs_dec_source_v4_last"] = source
        st.session_state["ngs_dec_method_v4_last"] = level

        precheck = precheck_decode_payload(
            payload,
            scheme_name=scheme_name,
            mode_codec=mode_codec,
            seed="rn",
            init_dimer=init_dimer,
            whiten=False,
            remove_leading_one=True,
            target_gc=0.50,
            w_gc=0.0,
            w_motif=0.0,
            ks=(4, 6),
        )
        st.session_state["ngs_dec_precheck_v4"] = precheck

        if level == "Manifest-guided rescue" and precheck.get("ok"):
            final_stats, restored = decode_dna_payload_to_file(
                payload,
                scheme_name=scheme_name,
                mode_codec=mode_codec,
                seed="rn",
                init_dimer=init_dimer,
                whiten=False,
                remove_leading_one=True,
                target_gc=0.50,
                w_gc=0.0,
                w_motif=0.0,
                ks=(4, 6),
                preferred_stem=preferred_stem,
            )
            st.session_state["ngs_dec_final_v4"] = final_stats
            st.session_state["ngs_dec_restored_v4"] = restored
        else:
            st.session_state["ngs_dec_final_v4"] = None
            st.session_state["ngs_dec_restored_v4"] = None

    payload = st.session_state.get("ngs_dec_payload_v4", "")
    report = st.session_state.get("ngs_dec_report_v4")
    precheck = st.session_state.get("ngs_dec_precheck_v4")
    final_stats = st.session_state.get("ngs_dec_final_v4")
    restored = st.session_state.get("ngs_dec_restored_v4")

    if report:
        with st.container(border=True):
            title = "### 📈 Recovery Summary" if (report.get("mode") == "standard_test") else "### 🛟 Recovery Summary"
            st.markdown(title)
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.metric("Total Reads", f"{int(report['total_reads']):,}")
            with c2:
                st.metric("Expected Fragments", f"{int(report['total_fragments']):,}")
            with c3:
                st.metric("Recovered Fragments", f"{int(report['recovered_fragments']):,}")
            with c4:
                st.metric("Recovery", f"{float(report['recovery_fraction']) * 100:.1f}%")
            if report.get("mode") == "advanced_rescue":
                d1, d2, d3 = st.columns(3)
                with d1:
                    st.metric("Rescued Fragments", f"{int(report.get('rescued_fragments') or 0):,}")
                with d2:
                    st.metric("Checksum-style Pass", f"{int(report.get('checksum_pass_like_fragments') or 0):,}")
                with d3:
                    st.metric("Uncertain", f"{int(report.get('uncertain_fragments') or 0):,}")
            report_df = pd.DataFrame(report["rows"]).copy()
            if "fragment_index" in report_df.columns:
                report_df.insert(0, "Fragment", report_df["fragment_index"].map(lambda x: f"Fragment {int(x)}"))
            if "fragment_name" in report_df.columns:
                report_df = report_df.drop(columns=["fragment_name"])
            report_df = report_df.rename(columns={
                "fragment_index": "Number",
                "supporting_reads": "Supporting Reads",
                "difficulty": "Difficulty",
                "score": "Score",
                "fragment_length_nt": "Fragment Length (nt)",
                "payload_length_nt": "Payload Length (nt)",
                "recovered": "Recovered",
                "rescued": "Rescued",
                "status": "Status",
                "mismatch_fraction_vs_manifest": "Mismatch Fraction",
            })
            st.dataframe(report_df, width="stretch", hide_index=True)

    if payload:
        with st.container(border=True):
            st.markdown("### 🧬 DNA Payload")
            _stats_block(payload)
            st.text_area("DNA payload preview", value=payload[:12000], height=220, key="ngs_dec_payload_preview_v4")
            st.download_button("Download DNA payload", data=(payload + "\n").encode("utf-8"), file_name="ngs_decode_payload.txt", width="stretch")

    if precheck:
        with st.container(border=True):
            st.markdown("### 🧪 Decode Check")
            ok = bool(precheck.get("ok"))
            if ok:
                st.success("The reconstructed DNA passes the decode and zlib integrity check.")
            else:
                st.error("The reconstructed DNA does not pass the decode and zlib integrity check.")
            st.json(precheck, expanded=False)

    if final_stats:
        with st.container(border=True):
            st.markdown("### 📥 Final Decode Result")
            d1, d2, d3 = st.columns(3)
            with d1:
                st.metric("Decoded bits", f"{int(final_stats.get('decoded_bits_len') or 0):,}")
            with d2:
                st.metric("Inner bytes", f"{int(final_stats.get('inner_bytes_len') or 0):,}")
            with d3:
                z = (final_stats.get("zlib") or {})
                st.metric("Zlib integrity", "OK" if z.get("integrity_ok") else "Failed")
            st.json(final_stats, expanded=False)
            if restored and os.path.exists(restored):
                st.success(f"Recovered file: {os.path.basename(restored)}")
                with open(restored, "rb") as f:
                    st.download_button("Download recovered file", data=f.read(), file_name=os.path.basename(restored), width="stretch")
