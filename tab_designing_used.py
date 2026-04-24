# tab_designing.py (Universal Step-1 tab)
# Streamlit UI port of "Universal DNA Compressor + Sequencing-Aware DNA Storage System"
# Focus: Step-1 (file <-> single DNA string, headerless at medium level).
#
# Backend:
#   - dna_codec.py (UNCHANGED)
#   - pipelines_v2.py, compressors_v2.py, utils_bits_v2.py (Step-1 v2 implementation)
#
# Run:
#   pip install -r requirements.txt
#   streamlit run tab_designing_universal_streamlit.py
#
# Notes (Windows):
#   - For Mode 3 audio/video conversions, ffmpeg must be installed and available on PATH
#     OR disable "allow_ffmpeg" to force keep-bytes.
#
import os
import re
import io
import json
import uuid
import tempfile
import subprocess
import shutil
import difflib
import hashlib
from typing import Any, Dict, Optional, Tuple, List
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import math
from PIL import Image
import dna_codec
from compressors_v2 import restore_rep, detect_domain
from pipelines_v2 import run_mode0_raw, run_mode1_zip, run_mode2_zip_store, run_mode3_domain, run_mode3_best
from utils_bits_v2 import (
    detect_magic,
    ensure_dir,
    safe_basename,
    zlib_inflate_until_eof,
    bitstring_to_bytes,
)
# ----------------------------
# Global dirs
# ----------------------------
ensure_dir("jobs")
ensure_dir("recovery_out")
ensure_dir(".streamlit_tmp_uploads")
# Refined UI styling
st.markdown("""
<style>
.stTabs [data-baseweb="tab"] p { font-size: 24px !important; font-weight: 700 !important; }
[data-testid="stHorizontalBlock"] .stTabs [data-baseweb="tab"] p { font-size: 18px !important; font-weight: 600 !important; }
:root {
    --dna-card-bg: linear-gradient(180deg, rgba(250,252,255,0.98), rgba(244,248,252,0.98));
    --dna-card-border: rgba(90, 125, 156, 0.22);
}
div[data-testid="stVerticalBlockBorderWrapper"] {
    background: var(--dna-card-bg);
    border: 1px solid var(--dna-card-border) !important;
    border-radius: 16px !important;
    box-shadow: 0 8px 24px rgba(15, 23, 42, 0.045);
}
div[data-testid="stMetric"] {
    background: rgba(255,255,255,0.72);
    border: 1px solid rgba(90,125,156,0.14);
    border-radius: 12px;
    padding: 0.2rem 0.55rem 0.15rem 0.55rem;
}
h3, h4 {
    margin-bottom: 0.2rem !important;
}
.block-container {
    padding-top: 1.15rem;
    padding-bottom: 1rem;
}
</style>
""", unsafe_allow_html=True)
# ----------------------------
# Helpers
# ----------------------------
def _compute_psnr_ssim(img_a_path: str, img_b_path: str):
    """Compute PSNR and (optionally) SSIM for two images. Returns (psnr, ssim)."""
    try:
        a = Image.open(img_a_path).convert("RGB")
        b = Image.open(img_b_path).convert("RGB")
        if a.size != b.size:
            b = b.resize(a.size)
        a_np = np.asarray(a).astype(np.float32)
        b_np = np.asarray(b).astype(np.float32)
        mse = float(np.mean((a_np - b_np) ** 2))
        if mse <= 1e-12:
            psnr = 99.0
        else:
            psnr = 20.0 * math.log10(255.0 / math.sqrt(mse))
        ssim_val = None
        try:
            from skimage.metrics import structural_similarity as ssim
            ssim_c = []
            for ch in range(3):
                ssim_c.append(float(ssim(a_np[:, :, ch], b_np[:, :, ch], data_range=255.0)))
            ssim_val = float(np.mean(ssim_c))
        except Exception:
            ssim_val = None
        return float(psnr), (float(ssim_val) if ssim_val is not None else None)
    except Exception:
        return None, None
def _ui_json_expander(data: dict, title: str = "Details (JSON)", height: int = 260):
    """UI helper: show JSON in a collapsible panel with fixed-height scroll to avoid pushing other widgets."""
    try:
        s = json.dumps(data, indent=2, ensure_ascii=False)
    except Exception:
        s = str(data)
    with st.expander(title, expanded=False):
        st.text_area(" ", value=s, height=height)
def _fmt_bytes(n: Optional[int]) -> str:
    if n is None:
        return "—"
    x = float(n)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if x < 1024.0 or unit == "TB":
            return f"{x:.2f} {unit}" if unit != "B" else f"{int(x)} B"
        x /= 1024.0
    return f"{x:.2f} TB"
def _guess_domain_from_upload(uploaded_file) -> str:
    if uploaded_file is None:
        return "—"
    try:
        raw = uploaded_file.getvalue()
    except Exception:
        raw = b""
    name = (uploaded_file.name or "upload.bin")
    try:
        return detect_domain(name, raw)
    except Exception:
        return "other"
def _load_dna_from_report(report: Dict[str, Any]) -> str:
    try:
        dna_path = ((report.get("artifacts", {}) or {}).get("dna_txt"))
        if dna_path and os.path.exists(dna_path):
            return dna_codec.clean_dna_text(open(dna_path, "r", encoding="utf-8", errors="ignore").read())
    except Exception:
        pass
    return ""
def _homopolymer_profile(dna: str) -> Dict[str, Any]:
    dna = dna_codec.clean_dna_text(dna) if hasattr(dna_codec, "clean_dna_text") else "".join([c for c in (dna or "").upper() if c in "ACGT"])
    if not dna:
        return {
            "longest": 0,
            "homo_count": 0,
            "count_ge2": 0,
            "count_ge3": 0,
            "count_ge4": 0,
            "total_runs": 0,
            "exact_len_1": 0,
            "exact_len_2": 0,
            "exact_len_3": 0,
            "exact_len_4": 0,
            "exact_len_ge5": 0,
        }
    runs = []
    cur = 1
    for i in range(1, len(dna)):
        if dna[i] == dna[i - 1]:
            cur += 1
        else:
            runs.append(cur)
            cur = 1
    runs.append(cur)
    return {
        "longest": max(runs),
        "homo_count": sum(1 for r in runs if r >= 2),
        "count_ge2": sum(1 for r in runs if r >= 2),
        "count_ge3": sum(1 for r in runs if r >= 3),
        "count_ge4": sum(1 for r in runs if r >= 4),
        "total_runs": len(runs),
        "exact_len_1": sum(1 for r in runs if r == 1),
        "exact_len_2": sum(1 for r in runs if r == 2),
        "exact_len_3": sum(1 for r in runs if r == 3),
        "exact_len_4": sum(1 for r in runs if r == 4),
        "exact_len_ge5": sum(1 for r in runs if r >= 5),
    }
def _get_homopolymer_profile_from_report(report: Dict[str, Any]) -> Dict[str, Any]:
    dna_meta = (report.get("dna", {}) or {}) if isinstance(report, dict) else {}
    hp = (dna_meta.get("homopolymer") or {}) if isinstance(dna_meta, dict) else {}
    needed = ["longest", "homo_count", "count_ge2", "count_ge3", "count_ge4", "exact_len_1", "exact_len_2", "exact_len_3", "exact_len_4", "exact_len_ge5"]
    if all(k in hp for k in needed):
        return hp
    dna = _load_dna_from_report(report)
    calc = _homopolymer_profile(dna)
    merged = dict(calc)
    if isinstance(hp, dict):
        merged.update({k: v for k, v in hp.items() if v is not None})
    return merged
def _plot_homopolymer_profile(job: str, hp: Dict[str, Any]):
    labels = ["1-nt", "2-nt", "3-nt", "4-nt", ">=5-nt"]
    vals = [int(hp.get("exact_len_1") or 0), int(hp.get("exact_len_2") or 0), int(hp.get("exact_len_3") or 0), int(hp.get("exact_len_4") or 0), int(hp.get("exact_len_ge5") or 0)]
    fig = go.Figure()
    fig.add_trace(go.Bar(x=labels, y=vals, width=[0.5] * 5))
    fig.update_layout(title="Homopolymer statistics chart", height=300, margin=dict(l=20, r=20, t=60, b=40), showlegend=False)
    fig.update_xaxes(title="Homopolymer length")
    fig.update_yaxes(title="Count")
    return fig

def _plot_bits_grouped(job: str, input_bits: int, compressed_bits: int):
    # Two thin bars: Input vs Compressed
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=["Input", "Compressed"],
        y=[input_bits, compressed_bits],
        width=[0.35, 0.35],
    ))
    fig.update_layout(
        title="Bit-size (before DNA) — Input vs Compressed",
        barmode="group",
        bargap=0.55,
        height=320,
        margin=dict(l=20, r=20, t=60, b=40),
        showlegend=False,
    )
    fig.update_yaxes(title="bits")
    return fig
def _plot_compression_dial(job: str, input_bytes: int, output_bytes: int):
    if output_bytes <= 0 or input_bytes <= 0:
        factor = 0.0
    else:
        factor = float(input_bytes) / float(output_bytes)
    max_v = max(2.0, factor * 1.2, 10.0)
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=factor,
        number={"suffix": "×"},
        title={"text": "Input size / Output size"},
        gauge={"axis": {"range": [0, max_v]}, "bar": {"thickness": 0.25}},
    ))
    fig.update_layout(height=320, margin=dict(l=20, r=20, t=60, b=20))
    return fig

def _plot_psnr_bar(job: str, psnr: float):
    fig = go.Figure()
    fig.add_trace(go.Bar(x=[job], y=[psnr], width=[0.35], name="PSNR (dB)"))
    fig.add_hline(y=20.0, line_dash="dash", line_color="red", annotation_text="Practical alert line", annotation_position="top left")
    fig.update_layout(title="PSNR (dB)", height=280, margin=dict(l=20, r=20, t=60, b=40), showlegend=False)
    fig.update_yaxes(title="dB")
    return fig

def _plot_ssim_bar(job: str, ssim_val: float):
    fig = go.Figure()
    fig.add_trace(go.Bar(x=[job], y=[ssim_val], width=[0.35], name="SSIM"))
    fig.add_hline(y=0.80, line_dash="dash", line_color="red", annotation_text="Practical alert line", annotation_position="top left")
    fig.update_layout(title="SSIM", height=280, margin=dict(l=20, r=20, t=60, b=40), showlegend=False)
    fig.update_yaxes(title="SSIM", range=[0, 1])
    return fig

def _render_benchmark_analysis(report: Dict[str, Any]):
    bench = ((report or {}).get("benchmark", {}) or {}) if isinstance(report, dict) else {}
    analysis = ((report or {}).get("analysis", {}) or {}) if isinstance(report, dict) else {}
    rows = list(bench.get("top3") or [])
    if not rows:
        st.info("No benchmark rows available for this run.")
        return

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Domain", str(bench.get("domain") or ((report.get("rep", {}) or {}).get("meta", {}) or {}).get("domain") or "—"))
    with c2:
        st.metric("Selected method", str(analysis.get("selected_method") or bench.get("selected_method") or "—"))
    with c3:
        st.metric("Candidates tested", str(int(analysis.get("candidates_tested") or bench.get("candidates_tested") or 0)))
    with c4:
        st.metric("Benchmark metric", "Output/Input")

    chart_df = pd.DataFrame([
        {
            "Method": str(r.get("method") or r.get("policy") or "candidate"),
            "Benchmark": float(r.get("benchmark_percent") or 0.0),
        }
        for r in rows
    ])
    fig = go.Figure()
    fig.add_trace(go.Bar(x=chart_df["Method"], y=chart_df["Benchmark"], width=[0.45] * len(chart_df)))
    fig.update_layout(title="Top 3 Benchmark Candidates", height=300, margin=dict(l=20, r=20, t=55, b=40), showlegend=False)
    fig.update_yaxes(title="Output / Input (%)")
    st.plotly_chart(fig, width="stretch", key=f"benchmark_top3_{(report or {}).get('job_uuid','empty')}")

    table_rows = []
    for r in rows:
        table_rows.append({
            "Rank": int(r.get("rank") or 0),
            "Method": str(r.get("method") or r.get("policy") or "candidate"),
            "Policy": str(r.get("policy") or "—"),
            "Lossy": "Yes" if bool(r.get("lossy", False)) else "No",
            "Rep size": _fmt_bytes(int(r.get("rep_size_bytes") or 0)),
            "Output size": _fmt_bytes(int(r.get("zlib_size_bytes") or 0)),
            "Benchmark": (f"{float(r.get('benchmark_percent')):.2f}%" if r.get("benchmark_percent") is not None else "—"),
            "Saving": (f"{float(r.get('size_saving_pct')):.2f}%" if r.get("size_saving_pct") is not None else "—"),
            "Δ vs best": (f"{float(r.get('delta_vs_best_pct')):.2f}%" if r.get("delta_vs_best_pct") is not None else "—"),
            "Selected": "Yes" if bool(r.get("selected", False)) else "No",
        })
    st.dataframe(pd.DataFrame(table_rows), width="stretch", hide_index=True)

def _word_edit_distance(a_words, b_words) -> int:
    """Levenshtein distance at word level."""
    n, m = len(a_words), len(b_words)
    if n == 0:
        return m
    if m == 0:
        return n
    dp = list(range(m + 1))
    for i in range(1, n + 1):
        prev = dp[0]
        dp[0] = i
        for j in range(1, m + 1):
            cur = dp[j]
            cost = 0 if a_words[i - 1] == b_words[j - 1] else 1
            dp[j] = min(dp[j] + 1, dp[j - 1] + 1, prev + cost)
            prev = cur
    return dp[m]
def _is_text_like(path: str) -> bool:
    ext = (os.path.splitext(path)[1] or "").lower()
    if ext in {".txt",".md",".csv",".tsv",".json",".xml",".html",".htm",".py",".js",".css",".yaml",".yml",".log"}:
        return True
    try:
        data = open(path, "rb").read(4096)
        if b"\x00" in data:
            return False
        s = data.decode("utf-8", errors="replace")
        repl = s.count("\ufffd")
        return repl / max(1, len(s)) < 0.02
    except Exception:
        return False
def _read_text_words(path: str):
    try:
        s = open(path, "rb").read().decode("utf-8", errors="ignore")
    except Exception:
        return []
    return [w for w in re.split(r"\s+", s.strip()) if w]


def _read_text_content(path: str) -> str:
    try:
        return open(path, "rb").read().decode("utf-8", errors="ignore")
    except Exception:
        return ""

def _char_edit_distance(a: str, b: str, cap: int = 12000) -> int:
    a = (a or "")[:cap]
    b = (b or "")[:cap]
    n, m = len(a), len(b)
    if n == 0:
        return m
    if m == 0:
        return n
    dp = list(range(m + 1))
    for i in range(1, n + 1):
        prev = dp[0]
        dp[0] = i
        ai = a[i - 1]
        for j in range(1, m + 1):
            cur = dp[j]
            cost = 0 if ai == b[j - 1] else 1
            dp[j] = min(dp[j] + 1, dp[j - 1] + 1, prev + cost)
            prev = cur
    return dp[m]

def _sha256_file(path: str) -> Optional[str]:
    try:
        h = hashlib.sha256()
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b''):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None

def _ffprobe_json(path: str) -> Dict[str, Any]:
    try:
        cmd = [
            'ffprobe', '-v', 'error', '-print_format', 'json',
            '-show_streams', '-show_format', path
        ]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return json.loads(out.stdout)
    except Exception:
        return {}

def _safe_float(v: Any) -> Optional[float]:
    try:
        return float(v)
    except Exception:
        return None

def _parse_fraction(v: str) -> Optional[float]:
    try:
        if not v:
            return None
        if '/' in v:
            a, b = v.split('/', 1)
            b = float(b)
            if b == 0:
                return None
            return float(a) / b
        return float(v)
    except Exception:
        return None

def _audio_meta(path: str) -> Dict[str, Any]:
    info = _ffprobe_json(path)
    fmt = info.get('format', {}) or {}
    streams = info.get('streams', []) or []
    a = next((s for s in streams if s.get('codec_type') == 'audio'), {})
    return {
        'duration_sec': _safe_float(fmt.get('duration') or a.get('duration')),
        'sample_rate': int(a.get('sample_rate')) if a.get('sample_rate') else None,
        'channels': int(a.get('channels')) if a.get('channels') else None,
        'codec': a.get('codec_name'),
        'bit_rate': int(fmt.get('bit_rate')) if fmt.get('bit_rate') and str(fmt.get('bit_rate')).isdigit() else None,
        'container': fmt.get('format_name'),
    }

def _video_meta(path: str) -> Dict[str, Any]:
    info = _ffprobe_json(path)
    fmt = info.get('format', {}) or {}
    streams = info.get('streams', []) or []
    v = next((s for s in streams if s.get('codec_type') == 'video'), {})
    a = next((s for s in streams if s.get('codec_type') == 'audio'), {})
    return {
        'duration_sec': _safe_float(fmt.get('duration') or v.get('duration')),
        'width': int(v.get('width')) if v.get('width') else None,
        'height': int(v.get('height')) if v.get('height') else None,
        'fps': _parse_fraction(v.get('avg_frame_rate') or v.get('r_frame_rate')),
        'video_codec': v.get('codec_name'),
        'audio_codec': a.get('codec_name'),
        'bit_rate': int(fmt.get('bit_rate')) if fmt.get('bit_rate') and str(fmt.get('bit_rate')).isdigit() else None,
        'container': fmt.get('format_name'),
    }

def _make_temp_png(prefix: str) -> str:
    tmp_dir = os.path.join(tempfile.gettempdir(), 'dna_storage_analysis')
    os.makedirs(tmp_dir, exist_ok=True)
    return os.path.join(tmp_dir, f"{prefix}_{uuid.uuid4().hex}.png")

def _audio_spectrogram_image(path: str) -> Optional[str]:
    out = _make_temp_png('spectrogram')
    try:
        cmd = [
            'ffmpeg', '-y', '-v', 'error', '-i', path,
            '-lavfi', 'showspectrumpic=s=1400x320:legend=disabled:color=viridis',
            '-frames:v', '1', out
        ]
        subprocess.run(cmd, check=True, capture_output=True)
        return out if os.path.exists(out) else None
    except Exception:
        return None

def _extract_video_frame(path: str, time_sec: float) -> Optional[str]:
    out = _make_temp_png('frame')
    try:
        cmd = ['ffmpeg', '-y', '-v', 'error', '-ss', f'{max(0.0, float(time_sec)):.3f}', '-i', path, '-frames:v', '1', out]
        subprocess.run(cmd, check=True, capture_output=True)
        return out if os.path.exists(out) else None
    except Exception:
        return None

def _video_keyframe_metrics(in_path: str, out_path: str) -> Dict[str, Any]:
    meta_a = _video_meta(in_path)
    meta_b = _video_meta(out_path)
    da = meta_a.get('duration_sec') or 0.0
    db = meta_b.get('duration_sec') or 0.0
    dur = min(da, db)
    if dur <= 0:
        return {'psnr': None, 'ssim': None, 'frame_a': None, 'frame_b': None, 'sample_time_sec': None}
    sample_time = max(0.0, dur * 0.5)
    fa = _extract_video_frame(in_path, sample_time)
    fb = _extract_video_frame(out_path, sample_time)
    psnr, ssim = (None, None)
    if fa and fb:
        psnr, ssim = _compute_psnr_ssim(fa, fb)
    return {'psnr': psnr, 'ssim': ssim, 'frame_a': fa, 'frame_b': fb, 'sample_time_sec': sample_time}

def _image_diff_array(img_a_path: str, img_b_path: str):
    try:
        a = Image.open(img_a_path).convert('RGB')
        b = Image.open(img_b_path).convert('RGB')
        if a.size != b.size:
            b = b.resize(a.size)
        a_np = np.asarray(a).astype(np.int16)
        b_np = np.asarray(b).astype(np.int16)
        diff = np.abs(a_np - b_np).astype(np.uint8)
        return diff
    except Exception:
        return None

def _domain_for_report(report: Optional[Dict[str, Any]]) -> str:
    if not report:
        return '—'
    rep = (report.get('rep', {}) or {})
    rep_meta = (rep.get('meta', {}) or {}) if isinstance(rep, dict) else {}
    return rep_meta.get('domain') or rep_meta.get('detected_domain') or '—'

def _get_image_analysis_payload(report: Optional[Dict[str, Any]]):
    art = (report.get('artifacts', {}) or {}) if report else {}
    out = (report.get('output', {}) or {}) if report else {}
    in_path = art.get('input_original') if report else None
    out_path = out.get('restored_file') if report else None
    psnr = ssim_val = None
    diff = None
    dims_a = dims_b = mode_a = mode_b = None
    if report and isinstance(in_path, str) and isinstance(out_path, str) and os.path.exists(in_path) and os.path.exists(out_path):
        psnr, ssim_val = _compute_psnr_ssim(in_path, out_path)
        diff = _image_diff_array(in_path, out_path)
        try:
            ia = Image.open(in_path)
            ib = Image.open(out_path)
            dims_a = f"{ia.width} × {ia.height}"
            dims_b = f"{ib.width} × {ib.height}"
            mode_a = ia.mode
            mode_b = ib.mode
        except Exception:
            pass
    return {
        'in_path': in_path,
        'out_path': out_path,
        'psnr': psnr,
        'ssim': ssim_val,
        'diff': diff,
        'dims_a': dims_a,
        'dims_b': dims_b,
        'mode_a': mode_a,
        'mode_b': mode_b,
    }


def _render_image_validation(report: Optional[Dict[str, Any]], job: str):
    payload = _get_image_analysis_payload(report)
    psnr = payload['psnr']
    ssim_val = payload['ssim']
    with st.container(border=True):
        st.markdown('### Image Validation')
        r1, r2 = st.columns(2)
        with r1:
            st.metric('PSNR (dB)', f"{psnr:.3f}" if psnr is not None else '—')
            st.metric('Original size', payload['dims_a'] or '—')
            st.metric('Original mode', payload['mode_a'] or '—')
        with r2:
            st.metric('SSIM', f"{ssim_val:.5f}" if ssim_val is not None else '—')
            st.metric('Restored size', payload['dims_b'] or '—')
            st.metric('Restored mode', payload['mode_b'] or '—')
        c1, c2 = st.columns(2)
        with c1:
            st.plotly_chart(_plot_psnr_bar(job, float(psnr)) if psnr is not None else _plot_psnr_bar(job, 0.0), width='stretch', key=f'image_psnr_{job}')
        with c2:
            st.plotly_chart(_plot_ssim_bar(job, float(ssim_val)) if ssim_val is not None else _plot_ssim_bar(job, 0.0), width='stretch', key=f'image_ssim_{job}')


def _render_image_preview(report: Optional[Dict[str, Any]]):
    payload = _get_image_analysis_payload(report)
    with st.container(border=True):
        st.markdown('### Preview and Difference')
        c1, c2, c3 = st.columns(3)
        with c1:
            st.caption('Original')
            if payload['in_path'] and os.path.exists(payload['in_path']):
                st.image(payload['in_path'], width='stretch')
            else:
                _empty_card_spacer(160)
        with c2:
            st.caption('Restored')
            if payload['out_path'] and os.path.exists(payload['out_path']):
                st.image(payload['out_path'], width='stretch')
            else:
                _empty_card_spacer(160)
        with c3:
            st.caption('Absolute difference')
            if payload['diff'] is not None:
                st.image(payload['diff'], width='stretch', clamp=True)
            else:
                _empty_card_spacer(160)

def _render_text_analysis(report: Optional[Dict[str, Any]]):
    left, right = st.columns([0.95, 1.05], gap='small')
    art = (report.get('artifacts', {}) or {}) if report else {}
    out = (report.get('output', {}) or {}) if report else {}
    in_path = art.get('input_original') if report else None
    out_path = out.get('restored_file') if report else None
    a = _read_text_content(in_path) if report and in_path and os.path.exists(in_path) else ''
    b = _read_text_content(out_path) if report and out_path and os.path.exists(out_path) else ''
    exact = (a == b) if report else None
    cer = char_acc = None
    wer = word_acc = None
    char_ed = word_ed = None
    if report and a is not None and b is not None:
        char_ed = _char_edit_distance(a, b)
        ref_chars = max(1, len(a))
        cer = float(char_ed) / ref_chars
        char_acc = max(0.0, 1.0 - cer)
        aw = [w for w in re.split(r'\s+', a.strip()) if w]
        bw = [w for w in re.split(r'\s+', b.strip()) if w]
        word_ed = _word_edit_distance(aw, bw)
        ref_words = max(1, len(aw))
        wer = float(word_ed) / ref_words
        word_acc = max(0.0, 1.0 - wer)
    with left:
        with st.container(border=True):
            st.markdown('### Text Validation')
            c1, c2 = st.columns(2)
            with c1:
                st.metric('Exact match', _bool_text(exact))
                st.metric('Character accuracy', f"{char_acc:.4f}" if char_acc is not None else '—')
                st.metric('Word accuracy', f"{word_acc:.4f}" if word_acc is not None else '—')
            with c2:
                st.metric('Character error rate', f"{cer:.4f}" if cer is not None else '—')
                st.metric('Word error rate', f"{wer:.4f}" if wer is not None else '—')
                st.metric('Word edit distance', f"{int(word_ed):,}" if word_ed is not None else '—')
            st.metric('Character edit distance', f"{int(char_ed):,}" if char_ed is not None else '—')
    with right:
        with st.container(border=True):
            st.markdown('### Preview and Difference')
            c1, c2 = st.columns(2)
            with c1:
                st.caption('Original text')
                st.text_area('orig_text_preview', value=(a[:2500] if a else ''), height=220, label_visibility='collapsed', key=f'text_orig_{report.get("job_uuid") if report else "empty"}')
            with c2:
                st.caption('Restored text')
                st.text_area('rest_text_preview', value=(b[:2500] if b else ''), height=220, label_visibility='collapsed', key=f'text_rest_{report.get("job_uuid") if report else "empty"}')
            if a or b:
                diff_text = "\n".join(list(difflib.unified_diff((a[:1200]).splitlines(), (b[:1200]).splitlines(), lineterm=''))[:160])
                st.caption('Unified diff preview')
                st.code(diff_text or 'No visible differences in the sampled preview.', language='diff')
            else:
                _empty_card_spacer(140)

def _render_audio_analysis(report: Optional[Dict[str, Any]], job: str):
    left, right = st.columns([0.95, 1.05], gap='small')
    art = (report.get('artifacts', {}) or {}) if report else {}
    out = (report.get('output', {}) or {}) if report else {}
    in_path = art.get('input_original') if report else None
    out_path = out.get('restored_file') if report else None
    meta_a = _audio_meta(in_path) if report and in_path and os.path.exists(in_path) else {}
    meta_b = _audio_meta(out_path) if report and out_path and os.path.exists(out_path) else {}
    dur_diff = None
    if meta_a.get('duration_sec') is not None and meta_b.get('duration_sec') is not None:
        dur_diff = abs(float(meta_a['duration_sec']) - float(meta_b['duration_sec']))
    spec_a = _audio_spectrogram_image(in_path) if report and in_path and os.path.exists(in_path) else None
    spec_b = _audio_spectrogram_image(out_path) if report and out_path and os.path.exists(out_path) else None
    with left:
        with st.container(border=True):
            st.markdown('### Audio Validation')
            c1, c2 = st.columns(2)
            with c1:
                st.metric('Original duration (s)', f"{meta_a.get('duration_sec'):.2f}" if meta_a.get('duration_sec') is not None else '—')
                st.metric('Original sample rate', f"{int(meta_a.get('sample_rate')):,} Hz" if meta_a.get('sample_rate') else '—')
                st.metric('Original channels', str(meta_a.get('channels')) if meta_a.get('channels') else '—')
            with c2:
                st.metric('Restored duration (s)', f"{meta_b.get('duration_sec'):.2f}" if meta_b.get('duration_sec') is not None else '—')
                st.metric('Restored sample rate', f"{int(meta_b.get('sample_rate')):,} Hz" if meta_b.get('sample_rate') else '—')
                st.metric('Restored channels', str(meta_b.get('channels')) if meta_b.get('channels') else '—')
            s1, s2, s3 = st.columns(3)
            with s1:
                st.metric('Duration difference', f"{dur_diff:.3f} s" if dur_diff is not None else '—')
            with s2:
                st.metric('Sample rate match', _bool_text(meta_a.get('sample_rate') == meta_b.get('sample_rate') if meta_a and meta_b else None))
            with s3:
                st.metric('Channel match', _bool_text(meta_a.get('channels') == meta_b.get('channels') if meta_a and meta_b else None))
            if in_path and os.path.exists(in_path):
                st.audio(in_path)
            if out_path and os.path.exists(out_path):
                st.audio(out_path)
    with right:
        with st.container(border=True):
            st.markdown('### Preview')
            c1, c2 = st.columns(2)
            with c1:
                st.caption('Original audio')
                if spec_a and os.path.exists(spec_a):
                    st.image(spec_a, width='stretch')
                else:
                    _empty_card_spacer(180)
            with c2:
                st.caption('Restored audio')
                if spec_b and os.path.exists(spec_b):
                    st.image(spec_b, width='stretch')
                else:
                    _empty_card_spacer(180)

def _render_video_analysis(report: Optional[Dict[str, Any]], job: str):
    left, right = st.columns([0.95, 1.05], gap='small')
    art = (report.get('artifacts', {}) or {}) if report else {}
    out = (report.get('output', {}) or {}) if report else {}
    in_path = art.get('input_original') if report else None
    out_path = out.get('restored_file') if report else None
    meta_a = _video_meta(in_path) if report and in_path and os.path.exists(in_path) else {}
    meta_b = _video_meta(out_path) if report and out_path and os.path.exists(out_path) else {}
    dur_diff = None
    if meta_a.get('duration_sec') is not None and meta_b.get('duration_sec') is not None:
        dur_diff = abs(float(meta_a['duration_sec']) - float(meta_b['duration_sec']))
    frame_metrics = _video_keyframe_metrics(in_path, out_path) if report and in_path and out_path and os.path.exists(in_path) and os.path.exists(out_path) else {}
    with left:
        with st.container(border=True):
            st.markdown('### Video Validation')
            c1, c2 = st.columns(2)
            with c1:
                st.metric('Original duration (s)', f"{meta_a.get('duration_sec'):.2f}" if meta_a.get('duration_sec') is not None else '—')
                st.metric('Original resolution', f"{meta_a.get('width')} × {meta_a.get('height')}" if meta_a.get('width') and meta_a.get('height') else '—')
                st.metric('Original FPS', f"{meta_a.get('fps'):.3f}" if meta_a.get('fps') is not None else '—')
            with c2:
                st.metric('Restored duration (s)', f"{meta_b.get('duration_sec'):.2f}" if meta_b.get('duration_sec') is not None else '—')
                st.metric('Restored resolution', f"{meta_b.get('width')} × {meta_b.get('height')}" if meta_b.get('width') and meta_b.get('height') else '—')
                st.metric('Restored FPS', f"{meta_b.get('fps'):.3f}" if meta_b.get('fps') is not None else '—')
            m1, m2, m3 = st.columns(3)
            with m1:
                st.metric('Duration difference', f"{dur_diff:.3f} s" if dur_diff is not None else '—')
            with m2:
                st.metric('Keyframe PSNR (dB)', f"{frame_metrics.get('psnr'):.3f}" if frame_metrics.get('psnr') is not None else '—')
            with m3:
                st.metric('Keyframe SSIM', f"{frame_metrics.get('ssim'):.5f}" if frame_metrics.get('ssim') is not None else '—')
            if in_path and os.path.exists(in_path):
                st.video(in_path)
            if out_path and os.path.exists(out_path):
                st.video(out_path)
    with right:
        with st.container(border=True):
            st.markdown('### Preview')
            if frame_metrics.get('sample_time_sec') is not None:
                st.caption(f"Representative frame sampled at {frame_metrics.get('sample_time_sec'):.2f} s")
            c1, c2 = st.columns(2)
            with c1:
                st.caption('Original frame')
                if frame_metrics.get('frame_a') and os.path.exists(frame_metrics['frame_a']):
                    st.image(frame_metrics['frame_a'], width='stretch')
                else:
                    _empty_card_spacer(180)
            with c2:
                st.caption('Restored frame')
                if frame_metrics.get('frame_b') and os.path.exists(frame_metrics['frame_b']):
                    st.image(frame_metrics['frame_b'], width='stretch')
                else:
                    _empty_card_spacer(180)

def _render_other_analysis(report: Optional[Dict[str, Any]]):
    left, right = st.columns([0.95, 1.05], gap='small')
    art = (report.get('artifacts', {}) or {}) if report else {}
    out = (report.get('output', {}) or {}) if report else {}
    in_path = art.get('input_original') if report else None
    out_path = out.get('restored_file') if report else None
    sha_in = _sha256_file(in_path) if report and in_path and os.path.exists(in_path) else None
    sha_out = _sha256_file(out_path) if report and out_path and os.path.exists(out_path) else None
    with left:
        with st.container(border=True):
            st.markdown('### File Validation')
            st.metric('Exact byte match', _bool_text((sha_in == sha_out) if sha_in and sha_out else None))
            st.metric('Original file size', _fmt_bytes(os.path.getsize(in_path)) if in_path and os.path.exists(in_path) else '—')
            st.metric('Restored file size', _fmt_bytes(os.path.getsize(out_path)) if out_path and os.path.exists(out_path) else '—')
            st.write(f"**Original SHA-256:** `{sha_in or '—'}`")
            st.write(f"**Restored SHA-256:** `{sha_out or '—'}`")
    with right:
        with st.container(border=True):
            st.markdown('### File Summary')
            st.write(f"**Detected input type:** `{_value_text(_guess_domain_from_upload(type('obj', (), {'name': os.path.basename(in_path)})()) if in_path else None)}`")
            st.write(f"**Original path:** `{_value_text(in_path)}`")
            st.write(f"**Restored path:** `{_value_text(out_path)}`")
            if out_path and os.path.exists(out_path):
                st.download_button('Download restored file', data=open(out_path, 'rb').read(), file_name=os.path.basename(out_path), width='stretch', key=f'other_dl_{report.get("job_uuid") if report else "empty"}')
            else:
                _empty_card_spacer(180)

def _render_domain_analysis(report: Optional[Dict[str, Any]], job: str):
    domain = _domain_for_report(report)
    if domain == 'image':
        _render_image_analysis(report, job)
    elif domain == 'text':
        _render_text_analysis(report)
    elif domain == 'audio':
        _render_audio_analysis(report, job)
    elif domain == 'video':
        _render_video_analysis(report, job)
    else:
        _render_other_analysis(report)
def _empty_card_spacer(height: int = 220):
    st.markdown(f"<div style='height:{height}px'></div>", unsafe_allow_html=True)
def _bool_text(v) -> str:
    if v is None:
        return "—"
    return "True" if bool(v) else "False"
def _value_text(v, formatter=None):
    if v is None or v == "":
        return "—"
    if formatter is None:
        return str(v)
    try:
        return formatter(v)
    except Exception:
        return str(v)
def _render_summary_rows(rows):
    for label, value in rows:
        st.write(f"**{label}:** `{value}`")
def _project_label(report: Optional[Dict[str, Any]]) -> str:
    if not report:
        return "Project —"
    idx = report.get("project_index")
    return f"Project {idx}" if idx is not None else f"Project {str(report.get('job_uuid', ''))[:8]}"

def _friendly_data_type(report: Optional[Dict[str, Any]], detected_domain: str = "—") -> str:
    if not report:
        return detected_domain
    rep = (report.get("rep", {}) or {})
    rep_meta = (rep.get("meta", {}) or {}) if isinstance(rep, dict) else {}
    return rep_meta.get("domain") or rep_meta.get("detected_domain") or detected_domain

def _friendly_option_name(report: Optional[Dict[str, Any]]) -> str:
    if not report:
        return "—"
    rep = (report.get("rep", {}) or {})
    rep_meta = (rep.get("meta", {}) or {}) if isinstance(rep, dict) else {}
    return report.get("ui_option_name") or rep_meta.get("chosen_candidate") or rep_meta.get("policy") or report.get("mode") or "—"

def _friendly_mode_label(report: Optional[Dict[str, Any]]) -> str:
    if not report:
        return "—"
    branch = report.get("ui_branch_label") or ("Non-compression" if report.get("mode") == "mode0_raw" else "Compression")
    return f"{branch}: {_friendly_option_name(report)}"

def _status_text(report: Optional[Dict[str, Any]]) -> str:
    if not report:
        return "—"
    return "Complete" if str(report.get("status")) == "ok" else _value_text(report.get("status"))

def _all_init_dimers() -> List[str]:
    return [a + b for a in "ACGT" for b in "ACGT" if a != b]
def _render_ratio_card(report: Optional[Dict[str, Any]], chart_key: str):
    st.markdown("### Compression Ratio")
    inp = (report.get("input", {}) or {}) if report else {}
    zlb = (report.get("zlib_stream", {}) or {}) if report else {}
    input_bytes = int(inp.get("size_bytes") or 0)
    output_bytes = int(zlb.get("size_bytes") or 0)
    ratio = (output_bytes / input_bytes) if input_bytes else None
    c1, c2 = st.columns(2)
    with c1:
        st.metric("Input size", _fmt_bytes(input_bytes) if input_bytes else "—")
    with c2:
        st.metric("Output size / Input size", f"{ratio * 100:.2f}%" if ratio is not None else "—")
    if report and input_bytes and output_bytes:
        st.plotly_chart(_plot_compression_dial(str(report.get("job_uuid", ""))[:8], input_bytes, output_bytes), width="stretch", key=chart_key)
    else:
        _empty_card_spacer(250)

def _render_dna_card(report: Optional[Dict[str, Any]], dna_preview: str, preview_key: str):
    st.markdown("### DNA Length")
    dna_meta = (report.get("dna", {}) or {}) if report else {}
    hp = _get_homopolymer_profile_from_report(report) if report else _homopolymer_profile("")
    dna_len = int(dna_meta.get("dna_len_nt") or 0)
    gc = dna_meta.get("gc_fraction")
    st.caption(f"ATGC preview • Length: {dna_len:,} nt" if dna_len else "ATGC preview • Length: —")
    a1, a2, a3 = st.columns(3)
    with a1:
        st.metric("GC content", f"{float(gc):.4f}" if gc is not None else "—")
    with a2:
        st.metric("Longest homopolymer length", f"{int(hp.get('longest') or 0)}" if report else "—")
    with a3:
        st.metric("Homopolymer segments (2+ bases)", f"{int(hp.get('homo_count') or hp.get('count_ge2') or 0)}" if report else "—")
    st.text_area("DNA", value=dna_preview if report else "", height=260, label_visibility="collapsed", key=preview_key)

def _render_summary_card(report: Optional[Dict[str, Any]], detected_domain: str = "—"):
    st.markdown("### Summary")
    inp = (report.get("input", {}) or {}) if report else {}
    dna = (report.get("dna", {}) or {}) if report else {}
    rows = [("Mode", _friendly_mode_label(report)), ("Status", _status_text(report)), ("Data type", _friendly_data_type(report, detected_domain) if report else detected_domain), ("Input size", _value_text(int(inp.get("size_bytes") or 0) if report and inp.get("size_bytes") is not None else None, _fmt_bytes)), ("DNA sequence length", _value_text(int(dna.get("dna_len_nt") or 0) if report and dna.get("dna_len_nt") is not None else None, lambda v: f"{int(v):,} nt"))]
    _render_summary_rows(rows)
    with st.expander("Show encode_report.json", expanded=False):
        st.json(report if report else {})

def _render_download_card(report: Optional[Dict[str, Any]], prefix: str, reset_button_label: Optional[str] = None, reset_state_key: Optional[str] = None, uploader_reset: bool = False):
    st.markdown("### Download")
    if report:
        dna_path = (report.get("artifacts", {}) or {}).get("dna_txt")
        dna_bytes = open(dna_path, "rb").read() if dna_path and os.path.exists(dna_path) else b""
        report_bytes = json.dumps(report, indent=2, ensure_ascii=False).encode("utf-8")
    else:
        dna_bytes = b""
        report_bytes = b"{}\n"
    st.download_button("Download encode.dna", data=dna_bytes, file_name="encode.dna.txt", width="stretch", key=f"{prefix}_download_dna")
    st.download_button("Download encode_report.json", data=report_bytes, file_name="encode_report.json", width="stretch", key=f"{prefix}_download_report")
    if reset_button_label and reset_state_key:
        st.markdown("<div style='height:10px'></div>", unsafe_allow_html=True)
        if st.button(reset_button_label, width="stretch", key=f"{prefix}_reset"):
            st.session_state[reset_state_key] = False
            if uploader_reset:
                st.session_state.uploader_key += 1
            st.rerun()

def _render_input_card(uploaded_file, title: str = "Input File"):
    st.markdown(f"### {title}")
    detected_domain = _guess_domain_from_upload(uploaded_file)
    if uploaded_file is not None:
        size_bytes = len(uploaded_file.getvalue())
        c1, c2 = st.columns(2)
        with c1:
            st.metric("File name", uploaded_file.name)
        with c2:
            st.metric("Input size", _fmt_bytes(size_bytes))
        st.write(f"**Data type:** `{detected_domain}`")
    else:
        c1, c2 = st.columns(2)
        with c1:
            st.metric("File name", "—")
        with c2:
            st.metric("Input size", "—")
        st.write("**Data type:** `—`")
    return detected_domain

def _build_analysis_df(reports: list) -> "pd.DataFrame":
    rows = []
    for r in reports:
        inp = (r.get("input", {}) or {})
        rep = (r.get("rep", {}) or {})
        zlb = (r.get("zlib_stream", {}) or {})
        dna = (r.get("dna", {}) or {})
        out = (r.get("output", {}) or {})
        art = (r.get("artifacts", {}) or {})
        flags = (r.get("flags", {}) or {})
        rep_meta = ((rep.get("meta", {}) or {}) if isinstance(rep, dict) else {})
        input_bytes = int(inp.get("size_bytes") or 0)
        rep_bytes   = int(rep.get("size_bytes") or 0) if isinstance(rep, dict) else 0
        zlib_bytes  = int(zlb.get("size_bytes") or 0) if isinstance(zlb, dict) else 0
        input_bits = input_bytes * 8
        rep_bits   = rep_bytes * 8
        zlib_bits  = zlib_bytes * 8
        dna_nt = dna.get("dna_len_nt") if isinstance(dna, dict) else None
        bits_per_nt = dna.get("bits_per_nt_est") if isinstance(dna, dict) else None
        dna_bits_est = None
        try:
            if dna_nt is not None and bits_per_nt is not None:
                dna_bits_est = float(dna_nt) * float(bits_per_nt)
        except Exception:
            dna_bits_est = None
        rows.append({
            "job_uuid": r.get("job_uuid"),
            "mode": r.get("mode"),
            "status": r.get("status"),
            "domain": rep_meta.get("domain") or rep_meta.get("detected_domain"),
            "lossy": rep_meta.get("lossy"),
            "input_bits": input_bits,
            "rep_bits": rep_bits,
            "zlib_bits": zlib_bits,
            "dna_nt": dna_nt,
            "bits_per_nt": bits_per_nt,
            "dna_bits_est": dna_bits_est,
            "ratio_rep_over_input": (rep_bytes / input_bytes) if input_bytes else None,
            "ratio_zlib_over_input": (zlib_bytes / input_bytes) if input_bytes else None,
            "zlib_ok": (flags.get("zlib_integrity_ok") if isinstance(flags, dict) else None),
            "restored_file": out.get("restored_file") if isinstance(out, dict) else None,
            "input_original": art.get("input_original") if isinstance(art, dict) else None,
        })
    return pd.DataFrame(rows)
def _md_kv(d: Dict[str, Any], title: str = "Summary") -> str:
    lines = [f"### {title}"]
    for k, v in d.items():
        lines.append(f"- **{k}**: `{v}`")
    return "\n".join(lines)
def _render_step1_report(report: Dict[str, Any]) -> str:
    if not report:
        return "No report."
    mode = report.get("mode")
    status = report.get("status")
    err = report.get("error")
    inp = report.get("input", {})
    rep = report.get("rep", {})
    zlib_s = report.get("zlib_stream", {})
    dna = report.get("dna", {})
    infl = report.get("inflate", {})
    flags = report.get("flags", {})
    out = report.get("output", {})
    summary = {
        "mode": mode,
        "status": status,
        "error": err,
        "input_bytes": inp.get("size_bytes"),
        "rep_bytes": rep.get("size_bytes"),
        "zlib_stream_bytes": zlib_s.get("size_bytes"),
        "dna_len_nt": dna.get("dna_len_nt"),
        "bits_per_nt_est": round(dna.get("bits_per_nt_est") or 0, 6) if dna.get("bits_per_nt_est") is not None else None,
        "zlib_eof": infl.get("eof"),
        "unused_tail_bytes": infl.get("unused_tail_len_bytes"),
        "dna_roundtrip_ok": flags.get("dna_roundtrip_ok"),
        "zlib_integrity_ok": flags.get("zlib_integrity_ok"),
        "magic_detect_ok": flags.get("magic_detect_ok"),
        "restored_file": out.get("restored_file"),
    }
    return _md_kv(summary, title="Step-1 Run")
def _read_dna_text(dna_text: str, dna_file) -> str:
    s = ""
    if dna_text and dna_text.strip():
        s = dna_text.strip()
    elif dna_file is not None:
        # dna_file is UploadedFile
        s = dna_file.getvalue().decode("utf-8", errors="ignore")
    else:
        return ""
    if hasattr(dna_codec, "clean_dna_text"):
        return dna_codec.clean_dna_text(s)
    return "".join([c for c in s.upper() if c in "ACGT"])
def _save_uploaded_file(uploaded_file) -> str:
    """Save Streamlit UploadedFile to a temp path and return filepath."""
    ensure_dir(".streamlit_tmp_uploads")
    name = uploaded_file.name
    # sanitize name
    safe = "".join([c if c.isalnum() or c in "._-" else "_" for c in name])
    tmp_dir = os.path.join(".streamlit_tmp_uploads", str(uuid.uuid4()))
    ensure_dir(tmp_dir)
    path = os.path.join(tmp_dir, safe)
    with open(path, "wb") as f:
        f.write(uploaded_file.getvalue())
    return path
def _auto_pick_mode_for_file(path: str, raw: bytes) -> str:
    """
    Auto routes through the domain-aware pipeline for every supported domain,
    including document/archive/binary/other. Mode 2 stays as a legacy/internal fallback
    but is no longer the normal auto path.
    """
    try:
        dom = detect_domain(path, raw)
    except Exception:
        dom = "other"
    return "mode3_domain" if dom in {"image", "audio", "video", "text", "document", "archive", "binary", "other"} else "mode2_zip_store"
def _apply_ui_labels(report: Dict[str, Any], branch_label: str, option_label: str) -> Dict[str, Any]:
    history = st.session_state.setdefault("history", [])
    report["ui_branch_label"] = branch_label
    report["ui_option_name"] = option_label
    report["project_index"] = len(history)
    return report

def _run_step1_streamlit(
    uploaded_file,
    mode_choice: str,
    # dna codec params
    scheme_name: str,
    codec_mode: str,
    seed: str,
    init_dimer: str,
    prepend_one: bool,
    whiten: bool,
    target_gc: float,
    w_gc: float,
    w_motif: float,
    ks: Tuple[int,int],
    zlib_policy: str,
    # mode1 zip
    zip_level: int,
    # mode3 policies
    allow_ffmpeg: bool,
    image_policy: str,
    webp_quality: int,
    text_policy: str,
    audio_policy: str,
    opus_bitrate_kbps: int,
    video_policy: str,
    video_crf: int,
    quality_mode: str = "Lossy",
    benchmark_best: bool = True,
) -> Tuple[str, Dict[str, Any], str]:
    """
    Returns: (job_dir, report, dna_preview)
    """
    input_path = _save_uploaded_file(uploaded_file)
    raw = open(input_path, "rb").read()
    common_kwargs = dict(
        scheme_name=scheme_name,
        mode_codec=codec_mode,
        seed=seed,
        init_dimer=init_dimer,
        prepend_one=prepend_one,
        remove_leading_one=True,
        whiten=whiten,
        target_gc=float(target_gc),
        w_gc=float(w_gc),
        w_motif=float(w_motif),
        ks=(int(ks[0]), int(ks[1])),
        zlib_policy=zlib_policy,
    )
    # Resolve selection
    if mode_choice in {"Auto", "Auto Fixed"}:
        picked = _auto_pick_mode_for_file(input_path, raw)
        if picked == "mode3_domain":
            job_dir, report = run_mode3_domain(
                input_path,
                image_policy=image_policy,
                webp_quality=int(webp_quality),
                text_policy=text_policy,
                allow_external_ffmpeg=bool(allow_ffmpeg),
                audio_policy=audio_policy,
                opus_bitrate_kbps=int(opus_bitrate_kbps),
                video_policy=video_policy,
                video_crf=int(video_crf),
                **common_kwargs,
            )
        else:
            job_dir, report = run_mode2_zip_store(input_path, **common_kwargs)
    elif mode_choice in {"Auto Best", "Auto Benchmark"}:
        picked = _auto_pick_mode_for_file(input_path, raw)
        if picked == "mode3_domain":
            job_dir, report = run_mode3_best(
                input_path,
                image_policy=image_policy,
                webp_quality=int(webp_quality),
                text_policy=text_policy,
                allow_external_ffmpeg=bool(allow_ffmpeg),
                audio_policy=audio_policy,
                opus_bitrate_kbps=int(opus_bitrate_kbps),
                video_policy=video_policy,
                video_crf=int(video_crf),
                **common_kwargs,
                quality_mode=quality_mode,
            )
        else:
            job_dir, report = run_mode2_zip_store(input_path, **common_kwargs)
    elif mode_choice.startswith("Branch 1") or mode_choice == "mode0_raw":
        job_dir, report = run_mode0_raw(input_path, **common_kwargs)
    elif mode_choice.startswith("Mode 1") or mode_choice == "mode1_zip":
        job_dir, report = run_mode1_zip(input_path, zip_level=int(zip_level), **common_kwargs)
    elif mode_choice.startswith("Mode 2") or mode_choice == "mode2_zip_store":
        job_dir, report = run_mode2_zip_store(input_path, **common_kwargs)
    else:
        job_dir, report = run_mode3_domain(
            input_path,
            image_policy=image_policy,
            webp_quality=int(webp_quality),
            text_policy=text_policy,
            allow_external_ffmpeg=bool(allow_ffmpeg),
            audio_policy=audio_policy,
            opus_bitrate_kbps=int(opus_bitrate_kbps),
            video_policy=video_policy,
            video_crf=int(video_crf),
            **common_kwargs,
        )
    dna_preview = ""
    dna_path = report.get("artifacts", {}).get("dna_txt")
    if dna_path and os.path.exists(dna_path):
        dna_preview = open(dna_path, "r", encoding="utf-8", errors="ignore").read().strip()[:12000]
    return job_dir, report, dna_preview
def _decode_step1_streamlit(
    dna_text: str,
    dna_file,
    scheme_name: str,
    codec_mode: str,
    seed: str,
    init_dimer: str,
    whiten: bool,
    remove_leading_one: bool,
    target_gc: float,
    w_gc: float,
    w_motif: float,
    ks: Tuple[int,int],
    preferred_stem: str,
) -> Tuple[Dict[str, Any], Optional[str]]:
    dna = _read_dna_text(dna_text, dna_file)
    if not dna:
        return {"error": "Provide DNA text or dna.txt"}, None
    bits, _digits = dna_codec.decode_dna_to_bits(
        dna,
        scheme_name=scheme_name,
        mode=codec_mode,
        seed=seed,
        init_dimer=init_dimer,
        remove_leading_one=bool(remove_leading_one),
        whiten=bool(whiten),
        target_gc=float(target_gc),
        w_gc=float(w_gc),
        w_motif=float(w_motif),
        ks=(int(ks[0]), int(ks[1])),
    )
    decoded_buf, pad_bits = bitstring_to_bytes(bits, pad_to_byte=True)
    inner, infl = zlib_inflate_until_eof(decoded_buf)
    z_ok = bool(infl.get("eof")) and (infl.get("error") is None)
    out_dir = os.path.join("recovery_out", str(uuid.uuid4()))
    ensure_dir(out_dir)
    stem = safe_basename(preferred_stem or "restored", fallback="restored")
    restored_file, restore_meta = restore_rep(inner, out_dir=out_dir, preferred_stem=stem)
    mk = detect_magic(inner)
    stats = {
        "dna": {"nt": len(dna), "bits_per_nt_est": (len(bits) / len(dna)) if len(dna) else None},
        "decoded": {"bits_len": len(bits), "bytes_len": len(decoded_buf), "pad_bits_to_byte": pad_bits},
        "zlib": {
            "eof": bool(infl.get("eof")),
            "unused_tail_len_bytes": int(infl.get("unused_tail_len", 0) or 0),
            "error": infl.get("error"),
            "integrity_ok": z_ok,
        },
        "inner": {"bytes_len": len(inner), "magic": (mk.kind if mk else None)},
        "output": {"restored_file": restored_file, "restore_meta": restore_meta},
    }
    return stats, (restored_file if z_ok else None)
# ----------------------------
# UI
# ----------------------------
def render_designing():
    st.header("DNA Data Design")
    if "uploader_key" not in st.session_state:
        st.session_state.uploader_key = 0
    if "history" not in st.session_state:
        st.session_state.history = []

    tab_enc, tab_dec, tab_ana = st.tabs(["Encoding", "Decoding", "Analysis"])

    with tab_enc:
        sub_raw, sub_comp = st.tabs(["Non-compression", "Compression"])

        def codec_controls(prefix: str, default_codec_mode: str = "TABLE"):
            st.markdown("### Initial Parameter")
            mode_label = st.selectbox("Mode", ["Simple Mapping", "Rule-base"], index=(0 if default_codec_mode == "SIMPLE" else 1), key=f"{prefix}_codec_mode")
            codec_mode = "SIMPLE" if mode_label == "Simple Mapping" else "TABLE"
            scheme_name = "RINF_B16"
            init_dimer = "TA"
            if codec_mode == "TABLE":
                c1, c2 = st.columns(2)
                with c1:
                    scheme_name = st.selectbox("DNA rule", ["R0_B9", "R1_B12", "R2_B15", "RINF_B16"], index=3, key=f"{prefix}_scheme")
                with c2:
                    dimers = _all_init_dimers()
                    init_dimer = st.selectbox("Initial dimer", dimers, index=dimers.index("TA"), key=f"{prefix}_init")
            return {"scheme_name": scheme_name, "codec_mode": codec_mode, "seed": "rn", "init_dimer": init_dimer, "prepend_one": (codec_mode != "SIMPLE"), "whiten": False, "target_gc": 0.50, "w_gc": 0.0, "w_motif": 0.0, "ks": (4, 6), "zlib_policy": "auto", "mode_label": mode_label}

        with sub_raw:
            st.subheader("Non-compression")
            left_col, mid_col, right_col = st.columns([1.0, 1.05, 1.0], gap="small")
            with left_col:
                with st.container(border=True):
                    up = st.file_uploader("Upload file", type=None, key=f"raw_u_{st.session_state.uploader_key}")
                    raw_domain = _render_input_card(up, "Input File")
                st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
                with st.container(border=True):
                    params = codec_controls("raw", default_codec_mode="SIMPLE")
                    if st.button("Run Encoding", type="primary", width="stretch", disabled=(up is None), key="raw_run_btn"):
                        st.session_state["run_raw"] = True
            report = None
            dna_preview = ""
            if up and st.session_state.get("run_raw"):
                try:
                    _job_dir, report, dna_preview = _run_step1_streamlit(up, "mode0_raw", params["scheme_name"], params["codec_mode"], params["seed"], params["init_dimer"], params["prepend_one"], params["whiten"], params["target_gc"], params["w_gc"], params["w_motif"], params["ks"], params["zlib_policy"], zip_level=6, allow_ffmpeg=False, image_policy="keep", webp_quality=80, text_policy="keep", audio_policy="keep", opus_bitrate_kbps=64, video_policy="keep", video_crf=28)
                    report = _apply_ui_labels(report, "Non-compression", params["mode_label"])
                    st.session_state.history.append(report)
                    st.session_state["run_raw"] = False
                except Exception as e:
                    st.error(f"Error: {e}")
            with mid_col:
                with st.container(border=True):
                    _render_ratio_card(report, chart_key=f"raw_ratio_chart_{(report or {}).get('job_uuid', 'empty')}")
                st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
                with st.container(border=True):
                    _render_dna_card(report, dna_preview, preview_key=f"raw_dna_preview_{(report or {}).get('job_uuid', 'empty')}")
            with right_col:
                with st.container(border=True):
                    _render_summary_card(report, raw_domain)
                st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
                with st.container(border=True):
                    _render_download_card(report, prefix=f"raw_{(report or {}).get('job_uuid', 'empty')}", reset_button_label="Reset", reset_state_key="run_raw", uploader_reset=True)

        with sub_comp:
            st.subheader("Compression")
            left_col, mid_col, right_col = st.columns([1.0, 1.05, 1.0], gap="small")
            with left_col:
                with st.container(border=True):
                    up = st.file_uploader("Upload file", type=None, key=f"comp_u_{st.session_state.uploader_key}")
                    detected_domain = _render_input_card(up, "Input File")
                st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
                with st.container(border=True):
                    st.markdown("### Initial Parameter")
                    mode_choice = "Auto Benchmark"
                    st.info("Compression mode is fixed to Auto Benchmark.")
                    params = codec_controls("comp", default_codec_mode="TABLE")
                    if st.button("Run Compression + DNA Encoding", type="primary", width="stretch", disabled=(up is None), key="comp_run_btn"):
                        st.session_state["run_comp"] = True
            report = None
            dna_preview = ""
            if up and st.session_state.get("run_comp"):
                try:
                    _job_dir, report, dna_preview = _run_step1_streamlit(up, mode_choice, params["scheme_name"], params["codec_mode"], params["seed"], params["init_dimer"], params["prepend_one"], params["whiten"], params["target_gc"], params["w_gc"], params["w_motif"], params["ks"], params["zlib_policy"], zip_level=6, allow_ffmpeg=True, image_policy="webp_lossy", webp_quality=80, text_policy="gzip", audio_policy="opus_ogg", opus_bitrate_kbps=64, video_policy="mp4_h264", video_crf=28, quality_mode="Lossy", benchmark_best=(mode_choice == "Auto Benchmark"))
                    report = _apply_ui_labels(report, "Compression", _friendly_option_name(report))
                    st.session_state.history.append(report)
                    st.session_state["run_comp"] = False
                except Exception as e:
                    st.error(f"Error: {e}")
            with mid_col:
                with st.container(border=True):
                    _render_ratio_card(report, chart_key=f"comp_ratio_chart_{(report or {}).get('job_uuid', 'empty')}")
                st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
                with st.container(border=True):
                    _render_dna_card(report, dna_preview, preview_key=f"comp_dna_preview_{(report or {}).get('job_uuid', 'empty')}")
            with right_col:
                with st.container(border=True):
                    _render_summary_card(report, detected_domain)
                st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
                with st.container(border=True):
                    _render_download_card(report, prefix=f"comp_{(report or {}).get('job_uuid', 'empty')}", reset_button_label="Reset", reset_state_key="run_comp", uploader_reset=True)

    with tab_dec:
        st.subheader("Decoding")
        left, right = st.columns([1.0, 1.0], gap="small")
        with left:
            with st.container(border=True):
                st.markdown("### Recovered")
                dna_text = st.text_area("Paste DNA sequence", value="", height=180, key="dec_dna_text")
                dna_file = st.file_uploader("Or upload DNA text file", type=["txt", "dna", "fasta", "fa"], key="dec_dna_file")
                dec_params = codec_controls("dec", default_codec_mode="TABLE")
                preferred_stem = st.text_input("Output name", value="restored", key="dec_stem")
                if st.button("Run Decoding", type="primary", width="stretch", key="dec_run_btn"):
                    st.session_state["run_dec"] = True
                restored = None
                stats = None
                if st.session_state.get("run_dec"):
                    try:
                        stats, restored = _decode_step1_streamlit(dna_text, dna_file, dec_params["scheme_name"], dec_params["codec_mode"], dec_params["seed"], dec_params["init_dimer"], dec_params["whiten"], True, dec_params["target_gc"], dec_params["w_gc"], dec_params["w_motif"], dec_params["ks"], preferred_stem)
                        st.session_state["run_dec"] = False
                    except Exception as e:
                        st.error(f"Error: {e}")
                        stats, restored = {"error": str(e)}, None
                if restored and os.path.exists(restored):
                    ext = os.path.splitext(restored)[1].lower()
                    if ext in {".png", ".jpg", ".jpeg", ".webp", ".bmp", ".gif", ".tif", ".tiff"}:
                        st.image(restored, width="stretch")
                    elif ext in {".mp4", ".mov", ".avi", ".mkv", ".webm", ".m4v"}:
                        st.video(restored)
                    elif ext in {".wav", ".mp3", ".flac", ".ogg", ".opus", ".m4a", ".aac"}:
                        st.audio(restored)
                    elif _is_text_like(restored):
                        st.text_area("Recovered text", value=_read_text_content(restored), height=220, key="dec_text_preview")
                    st.download_button("Download recovered file", data=open(restored, "rb").read(), file_name=os.path.basename(restored), width="stretch", key=f"dec_download_{os.path.basename(restored)}")
                else:
                    _empty_card_spacer(260)
        with right:
            with st.container(border=True):
                st.markdown("### Summary")
                latest = st.session_state.history[-1] if st.session_state.history else None
                branch = latest.get("ui_branch_label") if latest else None
                option_name = _friendly_option_name(latest) if latest else "—"
                rows = [("Mode", f"{branch}: {option_name}" if branch else "—"), ("Status", "Complete" if restored else ((stats or {}).get("error") or "—")), ("Data type", _friendly_data_type(latest) if latest else (_guess_domain_from_upload(type("obj", (), {"name": os.path.basename(restored)})()) if restored else "—")), ("Input size", _fmt_bytes(int((latest or {}).get("input", {}).get("size_bytes") or 0)) if latest else (_fmt_bytes(os.path.getsize(restored)) if restored and os.path.exists(restored) else "—")), ("DNA sequence length", f"{int(((latest or {}).get('dna', {}) or {}).get('dna_len_nt') or ((stats or {}).get('dna', {}) or {}).get('nt') or 0):,} nt" if (latest or stats) else "—")]
                _render_summary_rows(rows)
                with st.expander("Show decode_report.json", expanded=False):
                    st.json(stats if stats else {})

    with tab_ana:
        st.subheader("Analysis")
        report = st.session_state.history[-1] if st.session_state.history else None
        st.markdown(f"**{_project_label(report)}**")
        r = report or {}
        inp = (r.get("input", {}) or {})
        zlb = (r.get("zlib_stream", {}) or {})
        dna_meta = (r.get("dna", {}) or {})
        hp = _get_homopolymer_profile_from_report(r) if report else _homopolymer_profile("")
        input_bytes = int(inp.get("size_bytes") or 0)
        output_bytes = int(zlb.get("size_bytes") or 0)
        ratio = (output_bytes / input_bytes) if input_bytes else None
        bits_per_nt = dna_meta.get("bits_per_nt_est")
        top1, top2, top3, top4 = st.columns(4)
        with top1:
            st.metric("Input size", _fmt_bytes(input_bytes) if input_bytes else "—")
        with top2:
            st.metric("Compressed size", _fmt_bytes(output_bytes) if output_bytes else "—")
        with top3:
            st.metric("Output / Input", f"{ratio * 100:.2f}%" if ratio is not None else "—")
        with top4:
            st.metric("DNA sequence length", f"{int(dna_meta.get('dna_len_nt') or 0):,} nt" if report else "—")
        middle_left, middle_right = st.columns([1.02, 0.98], gap="small")
        with middle_left:
            with st.container(border=True):
                st.markdown("### Compression Ratio")
                rr1, rr2, rr3 = st.columns(3)
                with rr1:
                    st.metric("Input size", _fmt_bytes(input_bytes) if input_bytes else "—")
                with rr2:
                    st.metric("Output size", _fmt_bytes(output_bytes) if output_bytes else "—")
                with rr3:
                    st.metric("Output / Input", f"{ratio * 100:.2f}%" if ratio is not None else "—")
                if report and input_bytes and output_bytes:
                    st.plotly_chart(_plot_compression_dial(str(report.get("job_uuid", ""))[:8], input_bytes, output_bytes), width="stretch", key=f"analysis_ratio_gauge_{report.get('job_uuid', 'empty')}")
                else:
                    _empty_card_spacer(250)
            if _domain_for_report(report) == "image":
                st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
                _render_image_preview(report)
        with middle_right:
            if _domain_for_report(report) == "image":
                _render_image_validation(report, _project_label(report))
            elif _domain_for_report(report) == "text":
                _render_text_analysis(report)
            elif _domain_for_report(report) == "audio":
                _render_audio_analysis(report, _project_label(report))
            elif _domain_for_report(report) == "video":
                _render_video_analysis(report, _project_label(report))
            else:
                _render_other_analysis(report)
        with st.container(border=True):
            st.markdown("### Benchmark Results")
            _render_benchmark_analysis(report)
        with st.container(border=True):
            st.markdown("### DNA Sequence Statistics")
            d1, d2, d3, d4 = st.columns(4)
            with d1:
                st.metric("GC content", f"{float(dna_meta.get('gc_fraction')):.4f}" if report and dna_meta.get("gc_fraction") is not None else "—")
            with d2:
                st.metric("Longest homopolymer length", f"{int(hp.get('longest') or 0)}" if report else "—")
            with d3:
                st.metric("Homopolymer segments (2+ bases)", f"{int(hp.get('homo_count') or hp.get('count_ge2') or 0)}" if report else "—")
            with d4:
                st.metric("Bits per nt", f"{float(bits_per_nt):.4f}" if bits_per_nt is not None else "—")
            st.markdown("**Homopolymers**")
            h1, h2, h3, h4, h5 = st.columns(5)
            with h1:
                st.metric("1-nt", f"{int(hp.get('exact_len_1') or 0)}" if report else "—")
            with h2:
                st.metric("2-nt", f"{int(hp.get('exact_len_2') or 0)}" if report else "—")
            with h3:
                st.metric("3-nt", f"{int(hp.get('exact_len_3') or 0)}" if report else "—")
            with h4:
                st.metric("4-nt", f"{int(hp.get('exact_len_4') or 0)}" if report else "—")
            with h5:
                st.metric(">=5-nt", f"{int(hp.get('exact_len_ge5') or 0)}" if report else "—")
            if report:
                st.plotly_chart(_plot_homopolymer_profile(_project_label(report), hp), width="stretch", key=f"ana_homo_{report.get('job_uuid', 'empty')}")
            else:
                _empty_card_spacer(240)
