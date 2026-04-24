
import streamlit as st
from tab_home import render_home
from tab_designing_ui_clean import render_designing
from tab_about import render_about
from tab_guide_v1 import render_guide
from tab_ngs_prep_pipeline_v3 import render_ngs_prep
from tab_sequencing_fastq_v2 import render_sequencing_fastq
from tab_ngs_decode_pipeline_v3 import render_ngs_decode
from tab_recovery_analysis_v5 import render_recovery_analysis

st.set_page_config(page_title="DNA Data Storage Tool", page_icon="DDSS_logo.png", layout="wide", initial_sidebar_state="expanded")
st.markdown("""
<style>
div.block-container {padding-top: 6.4rem; padding-bottom: 1.2rem; padding-left: 1.4rem; padding-right: 1.4rem; max-width: 1480px;}
header[data-testid="stHeader"] {height: 0.01rem;}
h1, h2, h3 {margin-top: 0.1rem; margin-bottom: 0.45rem;}
p {margin-bottom: 0.35rem;}
section[data-testid="stSidebar"] {border-right: 1px solid rgba(60,72,88,0.10);}
section[data-testid="stSidebar"] .block-container {padding-top: 0.9rem;}
div[data-testid="stVerticalBlockBorderWrapper"] {border-radius:14px; border:1px solid rgba(77,99,120,0.16); background:linear-gradient(180deg, rgba(248,250,252,0.98), rgba(242,246,250,0.98)); box-shadow:0 4px 14px rgba(15,23,42,0.04);}
div[data-testid="stVerticalBlockBorderWrapper"] > div {padding-top:0.5rem; padding-bottom:0.55rem;}
div[role="radiogroup"] {display:flex; flex-wrap:nowrap; gap:0.4rem; margin-top:2.8rem; margin-bottom:0.4rem;}
div[role="radiogroup"] > label {padding:8px 14px; border-radius:12px; border:1px solid rgba(77,99,120,0.18); font-size:15px; font-weight:700; background:rgba(243,247,250,0.96);}
div[role="radiogroup"] > label:has(input:checked) {border:2px solid rgba(44,123,229,0.95); background:rgba(44,123,229,0.10);}
</style>
""", unsafe_allow_html=True)
with st.sidebar:
    try: st.image("DDSS_logo.png", width=280)
    except Exception: pass
    st.title("DNA Storage Lab")
    st.info("System Version: 1.0.0-Stable")
    st.markdown("Use the tabs above to move step by step through the DNA storage pipeline.")
PAGES={
    "Home": render_home,
    "Encode": render_designing,
    "NGS Prep": render_ngs_prep,
    "Sequencing / FASTQ": render_sequencing_fastq,
    "NGS Decode": render_ngs_decode,
    "Recovery Analysis": render_recovery_analysis,
    "Guide": render_guide,
    "About": render_about,
}
if "main_page" not in st.session_state: st.session_state["main_page"]="Encode"
st.markdown("<div style='height:5.8rem;'></div>", unsafe_allow_html=True)
page=st.radio("Main navigation", options=list(PAGES.keys()), key="main_page", horizontal=True, label_visibility="collapsed")
PAGES[page]()
