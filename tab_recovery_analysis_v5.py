
from __future__ import annotations
import csv, difflib, hashlib, io, json, os
import pandas as pd, streamlit as st
from ngs_decode_v5 import dna_stats
from ngs_prep_v5 import clean_dna
from preview_utils_v1 import render_file_preview

def _seq_similarity(a:str,b:str):
    a=clean_dna(a); b=clean_dna(b)
    if not a or not b:
        return {"exact_match":a==b,"same_length":len(a)==len(b),"position_match_fraction":None,"sequence_similarity":None}
    n=min(len(a),len(b)); pos=sum(1 for i in range(n) if a[i]==b[i])
    return {"exact_match":a==b,"same_length":len(a)==len(b),"position_match_fraction":(pos/n if n else None),"sequence_similarity":difflib.SequenceMatcher(None,a,b).ratio(),"mismatch_count_over_overlap":n-pos}

def render_recovery_analysis():
    st.title("📊 Recovery Analysis")
    payload_before=clean_dna(st.session_state.get("pipeline_payload_before_fastq") or "")
    payload_after=clean_dna(st.session_state.get("pipeline_decode_payload") or "")
    report=st.session_state.get("pipeline_decode_report"); precheck=st.session_state.get("pipeline_decode_precheck"); final_stats=st.session_state.get("pipeline_decode_final"); restored=st.session_state.get("pipeline_decode_restored")
    source=st.session_state.get("pipeline_decode_source"); method=st.session_state.get("pipeline_decode_method"); protection=st.session_state.get("pipeline_decode_protection_mode") or st.session_state.get("pipeline_protection_mode")
    prep=st.session_state.get("pipeline_ngs_prep_result")
    if not any([payload_after, report, precheck, final_stats]):
        st.info("No recovery result is available yet. Run NGS Decode first, then return to this tab."); return
    with st.container(border=True):
        st.subheader("Recovery Summary")
        c1,c2,c3,c4,c5,c6=st.columns(6)
        with c1: st.metric("Recovery Method", method or "—")
        with c2: st.metric("Decode Source", source or "—")
        with c3: st.metric("Protection Mode", protection or "—")
        with c4: st.metric("Reads", f"{int(report.get('total_reads') or 0):,}" if report else "—")
        with c5: st.metric("Fragments Recovered", f"{int(report.get('recovered_fragments') or 0):,}" if report else "—")
        with c6:
            final_status="Recovered file available" if restored and os.path.exists(restored) else ("Decode-ready" if precheck and precheck.get("ok") else "Needs attention")
            st.metric("Final Status", final_status)
    with st.container(border=True):
        st.subheader("DNA Payload Comparison")
        left,right=st.columns(2)
        with left:
            st.markdown("**Original payload statistics**"); s=dna_stats(payload_before) if payload_before else {"length_nt":0,"gc_ratio":0.0,"max_homopolymer":0,"shannon_entropy":0.0}
            a1,a2,a3,a4=st.columns(4)
            with a1: st.metric("Length", f"{s['length_nt']:,} nt" if payload_before else "—")
            with a2: st.metric("GC Ratio", f"{s['gc_ratio']:.4f}" if payload_before else "—")
            with a3: st.metric("Max Homopolymer", f"{s['max_homopolymer']} nt" if payload_before else "—")
            with a4: st.metric("Entropy", f"{s['shannon_entropy']:.4f}" if payload_before else "—")
        with right:
            st.markdown("**Recovered payload statistics**"); s=dna_stats(payload_after) if payload_after else {"length_nt":0,"gc_ratio":0.0,"max_homopolymer":0,"shannon_entropy":0.0}
            b1,b2,b3,b4=st.columns(4)
            with b1: st.metric("Length", f"{s['length_nt']:,} nt" if payload_after else "—")
            with b2: st.metric("GC Ratio", f"{s['gc_ratio']:.4f}" if payload_after else "—")
            with b3: st.metric("Max Homopolymer", f"{s['max_homopolymer']} nt" if payload_after else "—")
            with b4: st.metric("Entropy", f"{s['shannon_entropy']:.4f}" if payload_after else "—")
        sim=_seq_similarity(payload_before,payload_after)
        s1,s2,s3,s4=st.columns(4)
        with s1: st.metric("Exact Match", "Yes" if sim["exact_match"] else "No")
        with s2: st.metric("Position Match", f"{(sim['position_match_fraction'] or 0)*100:.2f}%" if sim["position_match_fraction"] is not None else "—")
        with s3: st.metric("Sequence Similarity", f"{(sim['sequence_similarity'] or 0)*100:.2f}%" if sim["sequence_similarity"] is not None else "—")
        with s4: st.metric("Same Length", "Yes" if sim["same_length"] else "No")
    with st.container(border=True):
        st.subheader("Recovery vs Overhead")
        if prep is not None:
            summ=prep.summary or {}
            c1,c2,c3,c4=st.columns(4)
            with c1: st.metric("Total Library Length", f"{int(summ.get('total_library_nt') or 0):,} nt")
            with c2: st.metric("Expansion Factor", f"{float(summ.get('expansion_factor') or 0):.2f}×" if summ.get("expansion_factor") is not None else "—")
            with c3: st.metric("Parity Fragments", f"{int(summ.get('parity_fragments') or 0):,}")
            with c4: st.metric("Recovery Rate", f"{float(report.get('recovery_fraction') or 0)*100:.1f}%" if report else "—")
    with st.container(border=True):
        st.subheader("Final Decode Result")
        d1,d2,d3,d4=st.columns(4)
        with d1: st.metric("Decode Bits", f"{int(final_stats.get('decoded_bits_len') or 0):,}" if final_stats else "—")
        with d2:
            z=(final_stats.get("zlib") or {}) if isinstance(final_stats,dict) else {}
            st.metric("zlib Integrity", "OK" if z.get("integrity_ok") else ("Failed" if final_stats else "—"))
        with d3: st.metric("Restored File", os.path.basename(restored) if restored and os.path.exists(restored) else "—")
        with d4: st.metric("File Hash Match", "Unavailable")
        if final_stats: st.json(final_stats, expanded=False)
        if restored and os.path.exists(restored): render_file_preview(restored, "Recovered Data Preview")
    with st.container(border=True):
        st.subheader("Failure Diagnosis")
        items=[]
        if report:
            frac=float(report.get("recovery_fraction") or 0.0)
            if frac < 0.8: items.append("Sequencing or fragment recovery is the main bottleneck because many fragments were not recovered.")
            elif frac < 1.0: items.append("Some fragments were only partially recovered, so the payload should be interpreted cautiously.")
            else: items.append("All expected fragments were recovered.")
            if report.get("rescued_fragments"): items.append("Some fragments were reconstructed using protection logic after errors were detected.")
            rows=report.get("rows") or []
            crc_failed=sum(1 for r in rows if r.get("status")=="CRC Failed")
            if crc_failed: items.append(f"{crc_failed} fragment(s) failed the embedded CRC check.")
        if precheck:
            if precheck.get("ok"): items.append("The reconstructed DNA passes the decode and zlib precheck.")
            else: items.append("The reconstructed DNA does not pass the decode or zlib precheck, so the payload is still corrupted or incomplete.")
        if final_stats:
            z=final_stats.get("zlib") or {}
            if z.get("integrity_ok") and restored and os.path.exists(restored): items.append("The pipeline produced a restored file successfully.")
            elif not z.get("integrity_ok"): items.append("The pipeline reached the zlib stage, but integrity validation failed.")
        if not items: items.append("Run NGS Decode first to generate a recovery result.")
        for it in items: st.write(f"• {it}")
    with st.container(border=True):
        st.subheader("Export")
        payload={"recovery_summary":report,"decode_source":source,"recovery_method":method,"protection_mode":protection,"precheck":precheck,"final_decode":final_stats,"payload_similarity":_seq_similarity(payload_before,payload_after)}
        st.download_button("Download recovery report JSON", data=json.dumps(payload, indent=2, ensure_ascii=False).encode("utf-8"), file_name="recovery_report.json", width="stretch")
        buf=io.StringIO(); w=csv.writer(buf); w.writerow(["Metric","Original","Recovered"]); orig=dna_stats(payload_before) if payload_before else {}; rec=dna_stats(payload_after) if payload_after else {}
        for k,label in [("length_nt","Length (nt)"),("gc_ratio","GC Ratio"),("max_homopolymer","Max Homopolymer"),("shannon_entropy","Shannon Entropy")]: w.writerow([label, orig.get(k), rec.get(k)])
        st.download_button("Download comparison CSV", data=buf.getvalue().encode("utf-8"), file_name="recovery_comparison.csv", width="stretch")
        diff=difflib.unified_diff([payload_before[i:i+120]+"\n" for i in range(0,len(payload_before),120)],[payload_after[i:i+120]+"\n" for i in range(0,len(payload_after),120)],fromfile="original_payload",tofile="recovered_payload")
        diff_text="".join(diff) or "No payload diff is available.\n"
        st.download_button("Download payload diff text", data=diff_text.encode("utf-8"), file_name="payload_diff.txt", width="stretch")
        log=[]
        if report: log.append(f"Recovery report: {json.dumps(report, ensure_ascii=False)}")
        if precheck: log.append(f"Precheck: {json.dumps(precheck, ensure_ascii=False)}")
        if final_stats: log.append(f"Final decode: {json.dumps(final_stats, ensure_ascii=False)}")
        st.download_button("Download decode log", data=("\n".join(log)+"\n").encode("utf-8"), file_name="decode_log.txt", width="stretch")
