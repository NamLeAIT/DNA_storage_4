
from __future__ import annotations

import os
from pathlib import Path
import streamlit as st


def render_file_preview(path: str, title: str = "Recovered Data Preview"):
    st.markdown(f"### 👁️ {title}")
    if not path or not os.path.exists(path):
        st.info("No file is available for preview.")
        return

    ext = os.path.splitext(path)[1].lower()
    image_ext = {".png", ".jpg", ".jpeg", ".bmp", ".gif", ".webp", ".tif", ".tiff"}
    audio_ext = {".wav", ".mp3", ".ogg", ".opus", ".flac", ".m4a", ".aac"}
    video_ext = {".mp4", ".mov", ".mkv", ".avi", ".webm", ".m4v"}
    text_ext = {".txt", ".md", ".json", ".csv", ".tsv", ".log", ".xml", ".yaml", ".yml", ".html", ".htm", ".py", ".js", ".c", ".cpp", ".java"}
    pdf_ext = {".pdf"}

    st.caption(f"File: {os.path.basename(path)}")

    if ext in image_ext:
        st.image(path, use_container_width=True)
        return
    if ext in audio_ext:
        st.audio(path)
        return
    if ext in video_ext:
        st.video(path)
        return
    if ext in text_ext:
        try:
            text = Path(path).read_text(encoding="utf-8", errors="ignore")
        except Exception:
            text = ""
        if text:
            st.text_area("Text preview", value=text[:12000], height=280)
        else:
            st.info("This text-like file could not be previewed.")
        return
    if ext in pdf_ext:
        with open(path, "rb") as f:
            st.download_button("Download recovered PDF", data=f.read(), file_name=os.path.basename(path), width="stretch")
        st.info("Inline PDF preview is not shown here. Please download the recovered PDF to inspect it.")
        return

    try:
        size_b = os.path.getsize(path)
        st.info(f"Preview is not available for this file type. File size: {size_b:,} bytes.")
    except Exception:
        st.info("Preview is not available for this file type.")
