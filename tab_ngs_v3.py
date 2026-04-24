
from __future__ import annotations

import json
import os
import uuid
from dataclasses import asdict
from typing import Optional

import pandas as pd
import streamlit as st

from ngs_prep_v3 import (
    NGSPrepConfig,
    PrimerConfig,
    clean_dna,
    gc_ratio,
    max_homopolymer,
    prepare_ngs_library,
    shannon_entropy,
    export_ngs_library,
    reconstruct_payload_from_fragments,
)
from ngs_fastq_v2 import (
    FASTQSimulationConfig,
    export_simulated_fastq,
    simulate_fastq_from_library,
)
from utils_bits_v2 import ensure_dir


def _load_latest_dna_from_history() -> str:
    history = st.session_state.get("history", [])
    if not history:
        return ""
    latest = history[-1]
    dna_path = ((latest.get("artifacts", {}) or {}).get("dna_txt"))
    if dna_path and os.path.exists(dna_path):
        try:
            return clean_dna(open(dna_path, "r", encoding="utf-8", errors="ignore").read())
        except Exception:
            return ""
    return ""


def _load_uploaded_text(uploaded_file) -> str:
    if uploaded_file is None:
        return ""
    try:
        raw = uploaded_file.getvalue().decode("utf-8", errors="ignore")
        return clean_dna(raw)
    except Exception:
        return ""


def _stats_block(seq: str):
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Total Length", f"{len(seq):,} nt" if seq else "—")
    with c2:
        st.metric("GC Ratio", f"{gc_ratio(seq):.4f}" if seq else "—")
    with c3:
        st.metric("Max Homopolymer", f"{max_homopolymer(seq)} nt" if seq else "—")
    with c4:
        st.metric("Shannon Entropy", f"{shannon_entropy(seq):.4f}" if seq else "—")


def render_ngs():
    st.title("🧫 NGS Preparation Mode")
    st.caption("Prepare DNA sequences for Next-Generation Sequencing with fragment structure, scoring, export, and optional FASTQ simulation.")

    src_col, cfg_col = st.columns([1.0, 1.0], gap="small")

    with src_col:
        with st.container(border=True):
            st.markdown("### 📥 Input")
            source = st.radio(
                "DNA Sequence Source",
                ["From Latest Design Run", "Upload DNA File", "Enter Directly"],
                horizontal=True,
                label_visibility="visible",
            )

            seq = ""
            if source == "From Latest Design Run":
                seq = _load_latest_dna_from_history()
                if seq:
                    st.success(f"Loaded DNA sequence ({len(seq):,} nt) from latest Design run")
                else:
                    st.warning("No DNA sequence found in session history yet. Run Design first or use another source.")
            elif source == "Upload DNA File":
                up = st.file_uploader("Upload DNA text / FASTA file", type=["txt", "dna", "fasta", "fa"], key="ngs_upload")
                seq = _load_uploaded_text(up)
                if seq:
                    st.success(f"Loaded DNA sequence ({len(seq):,} nt) from uploaded file")
            else:
                direct = st.text_area("Enter DNA sequence directly", value="", height=180, key="ngs_direct_dna")
                seq = clean_dna(direct)
                if seq:
                    st.success(f"Loaded DNA sequence ({len(seq):,} nt) from direct input")

        st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)

        with st.container(border=True):
            st.markdown("### 📊 Input Sequence Statistics")
            _stats_block(seq)

    with cfg_col:
        with st.container(border=True):
            st.markdown("### ⚙️ Fragmentation Settings")

            c1, c2 = st.columns(2)
            with c1:
                payload_len = st.number_input("Payload Length (nt)", min_value=20, max_value=200, value=80, step=1)
                fragment_id_len = st.number_input("Fragment ID Length (nt)", min_value=0, max_value=24, value=8, step=1)
                checksum_len = st.number_input("Checksum Length (nt)", min_value=0, max_value=12, value=0, step=1)
            with c2:
                last_policy_label = st.selectbox(
                    "Last Fragment Handling",
                    ["Keep shorter", "Pad to full length"],
                    index=0,
                )
                sample_barcode = st.text_input("Sample Barcode (optional)", value="", help="Short barcode prepended before fragment ID.")
                target_gc = st.slider("Target GC", min_value=0.35, max_value=0.65, value=0.50, step=0.01)

            st.markdown("#### 🔬 Primer Sequences")
            primer_source = st.selectbox("Primer Source", ["Default (Illumina-compatible)", "Custom"], index=0)
            if primer_source == "Default (Illumina-compatible)":
                fwd = "ACACGACGCTCTTCCGATCT"
                rev = "AGATCGGAAGAGCACACGTCT"
            else:
                fwd = st.text_input("Forward Primer (5')", value="ACACGACGCTCTTCCGATCT")
                rev = st.text_input("Reverse Primer (3')", value="AGATCGGAAGAGCACACGTCT")

            if seq:
                full_fragments = len(seq) // int(payload_len)
                remainder = len(seq) % int(payload_len)
                total_fragments = full_fragments + (1 if remainder > 0 else 0)
                standard_len = len(clean_dna(fwd)) + len(clean_dna(sample_barcode)) + int(fragment_id_len) + int(payload_len) + int(checksum_len) + len(clean_dna(rev))
                last_payload_len = remainder if remainder > 0 else int(payload_len)
                if last_policy_label == "Pad to full length":
                    last_payload_len = int(payload_len)
                last_frag_len = len(clean_dna(fwd)) + len(clean_dna(sample_barcode)) + int(fragment_id_len) + int(last_payload_len) + int(checksum_len) + len(clean_dna(rev))

                st.markdown("#### 📐 Fragmentation Calculation")
                k1, k2, k3, k4 = st.columns(4)
                with k1:
                    st.metric("Full Fragments", f"{full_fragments:,}")
                with k2:
                    st.metric("Standard Fragment Length", f"{standard_len} nt")
                with k3:
                    st.metric("Remainder", f"{remainder} nt")
                with k4:
                    st.metric("Last Fragment Length", f"{last_frag_len} nt")
            else:
                st.info("Load a DNA sequence to preview fragmentation statistics.")

    if "ngs_result" not in st.session_state:
        st.session_state["ngs_result"] = None
    if "ngs_exports" not in st.session_state:
        st.session_state["ngs_exports"] = None
    if "ngs_fastq_exports" not in st.session_state:
        st.session_state["ngs_fastq_exports"] = None
    if "ngs_sim_result" not in st.session_state:
        st.session_state["ngs_sim_result"] = None
    if "ngs_latest_exports" not in st.session_state:
        st.session_state["ngs_latest_exports"] = None

    run_btn = st.button("Generate NGS Fragments", type="primary", disabled=(not seq), width="stretch")

    if run_btn and seq:
        cfg = NGSPrepConfig(
            payload_length_nt=int(payload_len),
            fragment_id_length_nt=int(fragment_id_len),
            sample_barcode=sample_barcode,
            checksum_length_nt=int(checksum_len),
            primer=PrimerConfig(forward=fwd, reverse=rev),
            last_fragment_policy=("pad_to_full" if last_policy_label == "Pad to full length" else "keep_shorter"),
            target_gc=float(target_gc),
        )
        result = prepare_ngs_library(seq, cfg)
        out_dir = os.path.join("jobs", "ngs", str(uuid.uuid4()))
        ensure_dir(out_dir)
        exports = export_ngs_library(result, out_dir=out_dir, prefix="ngs_library")

        st.session_state["ngs_result"] = result
        st.session_state["ngs_exports"] = exports
        st.session_state["ngs_fastq_exports"] = None
        st.session_state["ngs_sim_result"] = None
        st.session_state["ngs_latest_exports"] = {"prep_result": result, "prep_exports": exports, "payload_before_fastq": reconstruct_payload_from_fragments(result.fragments)}
        st.success(f"Generated {result.summary['total_fragments']:,} NGS fragments!")

    result = st.session_state.get("ngs_result")
    exports = st.session_state.get("ngs_exports")

    if result:
        st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)

        with st.container(border=True):
            st.markdown("### 📤 Output")
            st.markdown("#### 📊 Fragment Summary")
            s = result.summary
            easy = s["difficulty_distribution"]["easy"]
            med = s["difficulty_distribution"]["medium"]
            hard = s["difficulty_distribution"]["hard"]
            total = max(1, s["total_fragments"])

            o1, o2, o3 = st.columns(3)
            with o1:
                st.metric("Total Fragments", f"{s['total_fragments']:,}")
            with o2:
                st.metric("Fragment Length", f"{s['standard_fragment_length_nt']} nt")
            with o3:
                st.metric("Avg. Score", f"{s['avg_score']:.1f}/100")

            st.markdown("#### 🎯 Sequencing Difficulty Distribution")
            d1, d2, d3 = st.columns(3)
            with d1:
                st.success(f"Easy: {easy} ({easy/total*100:.1f}%)")
            with d2:
                st.warning(f"Medium: {med} ({med/total*100:.1f}%)")
            with d3:
                st.error(f"Hard: {hard} ({hard/total*100:.1f}%)")

            rows = [asdict(x) for x in result.fragments]
            df = pd.DataFrame(rows)

            table_cols = [
                "fragment_name", "fragment_index", "payload_length_nt", "fragment_length_nt",
                "gc_ratio", "max_homopolymer", "shannon_entropy", "score", "difficulty",
                "is_padded", "is_last_fragment"
            ]
            pretty_df = df[table_cols].copy()
            pretty_df["gc_ratio"] = pretty_df["gc_ratio"].map(lambda x: round(float(x), 4))
            pretty_df["shannon_entropy"] = pretty_df["shannon_entropy"].map(lambda x: round(float(x), 4))

            st.markdown("#### 📋 Fragment Table")
            st.dataframe(pretty_df, width="stretch", hide_index=True)

            st.markdown("#### 🔍 Fragment Detail View")
            selected_name = st.selectbox("Select fragment to view", options=[r.fragment_name for r in result.fragments], index=0)
            rec = next((x for x in result.fragments if x.fragment_name == selected_name), None)
            if rec is not None:
                st.markdown("**Sequence Structure**")
                st.code(
                    f"""5' ─┬─ Forward Primer ({len(rec.forward_primer)} nt)
    │   {rec.forward_primer}
    ├─ Sample Barcode ({len(rec.sample_barcode)} nt)
    │   {rec.sample_barcode or '—'}
    ├─ Fragment ID ({len(rec.fragment_id)} nt)
    │   {rec.fragment_id or '—'}
    ├─ Payload ({len(rec.payload)} nt)
    │   {rec.payload[:60]}{'...' if len(rec.payload) > 60 else ''}
    ├─ Checksum ({len(rec.checksum)} nt)
    │   {rec.checksum or '—'}
    └─ Reverse Primer ({len(rec.reverse_primer)} nt)
        {rec.reverse_primer}
─── 3'""",
                    language=None,
                )
                detail_df = pd.DataFrame([
                    {"Metric": "GC Ratio", "Value": f"{rec.gc_ratio*100:.1f}%", "Status": "✅" if 0.40 <= rec.gc_ratio <= 0.60 else "⚠️"},
                    {"Metric": "Max Homopolymer", "Value": f"{rec.max_homopolymer} nt", "Status": "✅" if rec.max_homopolymer <= 4 else "⚠️"},
                    {"Metric": "Shannon Entropy", "Value": f"{rec.shannon_entropy:.4f}", "Status": "✅" if rec.shannon_entropy >= 1.8 else "⚠️"},
                    {"Metric": "Score", "Value": f"{rec.score}/100", "Status": "✅" if rec.score >= 85 else ("⚠️" if rec.score >= 65 else "❌")},
                    {"Metric": "Difficulty", "Value": rec.difficulty, "Status": ""},
                    {"Metric": "Self-complementarity", "Value": str(rec.self_complementarity_score), "Status": "✅" if rec.self_complementarity_score <= 2 else "⚠️"},
                ])
                st.dataframe(detail_df, width="stretch", hide_index=True)
                st.markdown("**Full Sequence**")
                st.code(rec.full_sequence, language=None)

            if result.fragments and result.fragments[-1].is_last_fragment and result.fragments[-1].fragment_length_nt != result.fragments[0].fragment_length_nt:
                st.info(
                    f"Last fragment ({result.fragments[-1].fragment_name}) has different length: "
                    f"{result.fragments[-1].fragment_length_nt} nt (payload: {result.fragments[-1].payload_length_nt} nt)"
                )

        with st.container(border=True):
            st.markdown("### 💾 Download")
            if exports:
                if os.path.exists(exports["fasta"]):
                    st.download_button("Download FASTA", data=open(exports["fasta"], "rb").read(), file_name=os.path.basename(exports["fasta"]), width="stretch")
                if os.path.exists(exports["csv"]):
                    st.download_button("Download CSV Table", data=open(exports["csv"], "rb").read(), file_name=os.path.basename(exports["csv"]), width="stretch")
                if os.path.exists(exports["json"]):
                    st.download_button("Download Manifest JSON", data=open(exports["json"], "rb").read(), file_name=os.path.basename(exports["json"]), width="stretch")

        with st.container(border=True):
            st.markdown("### 🧪 Optional FASTQ Simulation")
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                coverage = st.number_input("Coverage", min_value=1, max_value=500, value=20, step=1)
            with c2:
                paired_end = st.checkbox("Paired-end", value=True)
            with c3:
                read_length = st.number_input("Read Length", min_value=20, max_value=300, value=150, step=1)
            with c4:
                mean_q = st.number_input("Mean Q", min_value=20, max_value=40, value=34, step=1)

            e1, e2, e3 = st.columns(3)
            with e1:
                sub_rate = st.number_input("Substitution Rate", min_value=0.0, max_value=0.05, value=0.003, step=0.001, format="%.3f")
            with e2:
                ins_rate = st.number_input("Insertion Rate", min_value=0.0, max_value=0.02, value=0.0002, step=0.0001, format="%.4f")
            with e3:
                del_rate = st.number_input("Deletion Rate", min_value=0.0, max_value=0.02, value=0.0002, step=0.0001, format="%.4f")

            sim_btn = st.button("Generate Simulated FASTQ", width="stretch")
            if sim_btn:
                sim_cfg = FASTQSimulationConfig(
                    coverage=int(coverage),
                    paired_end=bool(paired_end),
                    read_length_nt=int(read_length),
                    mean_q=int(mean_q),
                    sub_rate=float(sub_rate),
                    ins_rate=float(ins_rate),
                    del_rate=float(del_rate),
                )
                sim = simulate_fastq_from_library(result, sim_cfg)
                out_dir = os.path.join("jobs", "ngs", str(uuid.uuid4()))
                ensure_dir(out_dir)
                sim_exports = export_simulated_fastq(sim, out_dir=out_dir, prefix="simulated_reads")
                st.session_state["ngs_fastq_exports"] = sim_exports
                st.session_state["ngs_sim_result"] = sim
                st.session_state["ngs_latest_exports"] = {"prep_result": result, "prep_exports": exports, "sim_result": sim, "sim_exports": sim_exports, "payload_before_fastq": reconstruct_payload_from_fragments(result.fragments)}
                st.success("Simulated FASTQ generated successfully.")

            sim_exports = st.session_state.get("ngs_fastq_exports")
            if sim_exports:
                if "r1_fastq_gz" in sim_exports and os.path.exists(sim_exports["r1_fastq_gz"]):
                    st.download_button("Download R1 FASTQ.gz", data=open(sim_exports["r1_fastq_gz"], "rb").read(), file_name=os.path.basename(sim_exports["r1_fastq_gz"]), width="stretch")
                if "r2_fastq_gz" in sim_exports and os.path.exists(sim_exports["r2_fastq_gz"]):
                    st.download_button("Download R2 FASTQ.gz", data=open(sim_exports["r2_fastq_gz"], "rb").read(), file_name=os.path.basename(sim_exports["r2_fastq_gz"]), width="stretch")
                if "summary_json" in sim_exports and os.path.exists(sim_exports["summary_json"]):
                    st.download_button("Download FASTQ Summary JSON", data=open(sim_exports["summary_json"], "rb").read(), file_name=os.path.basename(sim_exports["summary_json"]), width="stretch")
