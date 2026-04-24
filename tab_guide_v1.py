
from __future__ import annotations

import streamlit as st


def render_guide():
    st.title("📚 Guide")
    st.caption("A practical guide to the DNA data storage workflow in this app.")

    top_tabs = st.tabs(["ℹ️ About", "📖 Instructions", "❓ Help"])

    with top_tabs[0]:
        with st.container(border=True):
            st.subheader("About the Platform")
            st.markdown(
                """
                **DNA Data Storage Platform v2** simulates an end-to-end DNA data storage workflow:

                1. **Encode** a source file into a DNA payload  
                2. **Prepare NGS fragments** with primers and fragment identifiers  
                3. **Simulate or import FASTQ reads**  
                4. **Recover DNA** from sequencing-style data  
                5. **Decode the recovered payload** back into a file
                """
            )

        with st.container(border=True):
            st.subheader("What this app is for")
            st.markdown(
                """
                - Explore how digital files can be stored as DNA sequences  
                - Compare compression choices before DNA encoding  
                - Inspect fragment libraries designed for sequencing  
                - Simulate sequencing errors and test recovery  
                - Study where recovery succeeds or fails in the pipeline
                """
            )

    with top_tabs[1]:
        mode = st.selectbox(
            "Select a workflow step",
            ["Encode", "NGS Preparation", "Sequencing and FASTQ", "NGS Decode"],
            index=0,
            key="guide_mode_select",
        )

        if mode == "Encode":
            with st.container(border=True):
                st.subheader("📤 Encode")
                st.markdown(
                    """
                    **Step 1 — Upload a file**  
                    Upload the original file you want to store in DNA.

                    **Step 2 — Review file information**  
                    Check the detected file type, file size, and any preview data shown by the app.

                    **Step 3 — Run the encoding pipeline**  
                    The app automatically benchmarks suitable compression methods, selects the best result, applies zlib framing, and encodes the result as a DNA payload.

                    **Step 4 — Review the output**  
                    Inspect DNA length, GC ratio, homopolymer statistics, and the selected compression method.

                    **Step 5 — Download the DNA payload**  
                    Save the DNA payload if you want to continue manually or use it outside the app.
                    """
                )

        elif mode == "NGS Preparation":
            with st.container(border=True):
                st.subheader("🧫 NGS Preparation")
                st.markdown(
                    """
                    **Step 1 — Choose a DNA input source**  
                    You can use the latest DNA payload from Encode, upload a DNA file, or paste a DNA sequence directly.

                    **Step 2 — Set fragment options**  
                    Choose payload length, fragment identifier length, checksum length, and last-fragment handling.

                    **Step 3 — Set primer sequences**  
                    Use the default Illumina-compatible primers or enter custom primers.

                    **Step 4 — Generate the fragment library**  
                    The app builds sequencing-ready fragments and scores them for GC balance, homopolymers, entropy, and structure.

                    **Step 5 — Export the library**  
                    Download the fragment FASTA file, fragment table, and manifest JSON for downstream sequencing or recovery.
                    """
                )

        elif mode == "Sequencing and FASTQ":
            with st.container(border=True):
                st.subheader("🧪 Sequencing and FASTQ")
                st.markdown(
                    """
                    **Simulation mode**  
                    Generate sequencing-style FASTQ reads from the fragment library.

                    **Inputs needed**  
                    - A fragment library produced by **NGS Preparation**

                    **Options**  
                    - Coverage  
                    - Single-end or paired-end reads  
                    - Read length  
                    - Mean Q-score  
                    - Substitution, insertion, and deletion rates

                    **Output files**  
                    - R1 FASTQ  
                    - R2 FASTQ (if paired-end is enabled)  
                    - FASTQ summary JSON
                    """
                )

        else:
            with st.container(border=True):
                st.subheader("📥 NGS Decode")
                st.markdown(
                    """
                    **Available input sources**

                    1. **DNA before sequencing**  
                       Use the clean DNA payload before FASTQ simulation or sequencing.

                    2. **Latest FASTQ run**  
                       Use the most recent FASTQ data generated inside the app.

                    3. **FASTQ files from disk**  
                       Use this when you want to test recovery from sequencing-style reads.

                       **Files needed**
                       - Manifest JSON from **NGS Preparation**
                       - R1 FASTQ
                       - Optional R2 FASTQ

                    4. **Recovered DNA text file**  
                       Use a DNA text file that has already been reconstructed outside the app.

                    **Recommended workflow**
                    - First run **Recover DNA from reads**
                    - Then run **Decode recovered DNA to file**
                    """
                )

    with top_tabs[2]:
        with st.container(border=True):
            st.subheader("Frequently needed files")
            st.markdown(
                """
                **Manifest JSON**  
                Exported by the **NGS Preparation** tab. It describes the fragment library and is required for recovery from FASTQ.

                **R1 FASTQ / R2 FASTQ**  
                Produced by the sequencing simulation step or imported from an external sequencing workflow.

                **DNA payload text**  
                Produced by the **Encode** tab or reconstructed during recovery.
                """
            )

        with st.container(border=True):
            st.subheader("Common questions")
            st.markdown(
                """
                **Why does recovery fail even when FASTQ exists?**  
                Recovery can fail if too many fragments are missing or corrupted, or if the reconstructed DNA does not pass zlib integrity checks.

                **Why do fragment scores change after NGS preparation?**  
                Fragment scores are computed on the full sequencing-ready fragment, not only on the clean payload. Primers, fragment identifiers, checksums, and padding can change GC balance and homopolymer behavior.

                **What should I do first if I am new to the app?**  
                Start with **Encode**, then **NGS Preparation**, then simulate FASTQ, and finally use **NGS Decode**.
                """
            )
