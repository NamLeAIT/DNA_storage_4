
from __future__ import annotations

import streamlit as st


def render_guide():
    st.title("📚 Guide")
    st.caption("A practical user guide for the sequencing-aware DNA data storage pipeline in this application.")

    tabs = st.tabs(["Overview", "Step-by-Step Workflow", "Files You Need", "Help and FAQs"])

    with tabs[0]:
        with st.container(border=True):
            st.subheader("What this platform does")
            st.markdown(
                """
                **DNA Data Storage Platform** is a practical workflow simulator for storing digital data in DNA and
                testing how well that data can be recovered after sequencing-style noise.

                The current version of the platform supports the following pipeline:

                1. **Encode** a source file into a DNA payload  
                2. **Prepare sequencing-ready fragments** in **NGS Prep**  
                3. **Simulate or import FASTQ reads** in **Sequencing / FASTQ**  
                4. **Recover DNA and decode it** in **NGS Decode**  
                5. **Assess recovery reliability** in **Recovery Analysis**
                """
            )

        with st.container(border=True):
            st.subheader("Main tabs and their roles")
            st.markdown(
                """
                - **Home** — quick overview of the platform and workflow  
                - **Encode** — converts a source file into a DNA payload using compression and DNA encoding  
                - **NGS Prep** — builds sequencing-ready fragment libraries with primers, fragment identifiers, and optional protection modes  
                - **Sequencing / FASTQ** — simulates sequencing-style FASTQ reads with optional errors, or loads external FASTQ files  
                - **NGS Decode** — reconstructs DNA from clean payloads, fragment libraries, or FASTQ data, and optionally restores the original file  
                - **Recovery Analysis** — evaluates whether recovery succeeded, how accurate it was, and where failures occurred  
                - **Guide** — instructions, file requirements, and usage tips  
                - **About** — platform description, project scope, and contact information
                """
            )

        with st.container(border=True):
            st.subheader("Protection modes")
            st.markdown(
                """
                In the current version, the platform supports multiple protection modes during **NGS Prep**:

                - **None** — no protection, used as a baseline  
                - **CRC** — fragment-level integrity checking for error detection  
                - **Parity Fragments** — simple redundancy across fragment groups for limited fragment recovery  
                - **Reed–Solomon** — stronger block-level protection for more robust recovery under sequencing errors

                Each run uses **one protection mode at a time**, so the recovery result can be interpreted clearly.
                """
            )

    with tabs[1]:
        step = st.selectbox(
            "Select a workflow step",
            [
                "Encode",
                "NGS Prep",
                "Sequencing / FASTQ",
                "NGS Decode",
                "Recovery Analysis",
            ],
            index=0,
            key="guide_step_select_v2",
        )

        if step == "Encode":
            with st.container(border=True):
                st.subheader("📤 Encode")
                st.markdown(
                    """
                    **Purpose**  
                    Convert an original source file into a DNA payload.

                    **How to use it**
                    1. Upload the original file  
                    2. Review the detected file information  
                    3. Run the encoding pipeline  
                    4. Inspect the selected compression result and DNA statistics  
                    5. Download the DNA payload if needed

                    **What happens in this step**
                    - The app benchmarks suitable compression methods automatically  
                    - The selected representation is wrapped using zlib framing  
                    - The result is converted into a DNA payload using the selected DNA mapping method

                    **Typical outputs**
                    - Selected compression method  
                    - DNA payload  
                    - DNA length  
                    - GC ratio  
                    - Homopolymer statistics
                    """
                )

        elif step == "NGS Prep":
            with st.container(border=True):
                st.subheader("🧫 NGS Prep")
                st.markdown(
                    """
                    **Purpose**  
                    Convert a clean DNA payload into a sequencing-ready fragment library.

                    **How to use it**
                    1. Choose the DNA payload source  
                    2. Set the fragment size and primer settings  
                    3. Select a protection mode  
                    4. Generate sequencing fragments  
                    5. Review fragment statistics and export the library

                    **What happens in this step**
                    - The DNA payload is split into fragments  
                    - Primers and fragment identifiers are added  
                    - The selected protection method is applied  
                    - The app evaluates fragment quality and storage overhead

                    **Typical outputs**
                    - Fragment FASTA  
                    - Manifest JSON  
                    - Fragment CSV table  
                    - Protection overhead and expansion statistics
                    """
                )

        elif step == "Sequencing / FASTQ":
            with st.container(border=True):
                st.subheader("🧪 Sequencing / FASTQ")
                st.markdown(
                    """
                    **Purpose**  
                    Generate sequencing-style FASTQ reads with optional noise, or load external FASTQ files.

                    **How to use it**
                    1. Choose whether to simulate sequencing or load external FASTQ  
                    2. If simulating, set coverage and error parameters  
                    3. Run the simulation  
                    4. Review the injected error summary and fragment-level error statistics  
                    5. Export the FASTQ files

                    **What happens in this step**
                    - The fragment library is converted into reads  
                    - Optional substitution, insertion, and deletion errors are injected  
                    - The output FASTQ files are labeled clearly as noisy sequencing data

                    **Typical outputs**
                    - R1 FASTQ  
                    - R2 FASTQ  
                    - FASTQ summary JSON  
                    - Fragment-level error summary
                    """
                )

        elif step == "NGS Decode":
            with st.container(border=True):
                st.subheader("📥 NGS Decode")
                st.markdown(
                    """
                    **Purpose**  
                    Reconstruct DNA from clean DNA, fragment libraries, or FASTQ reads, then optionally decode it back into a file.

                    **How to use it**
                    1. Choose the decode source  
                    2. Choose the recovery method  
                    3. Run recovery  
                    4. Review the recovered DNA payload  
                    5. Run final decode to restore the file if the payload passes integrity checks

                    **Important note**  
                    The app uses the stored or manifest-based parameters automatically where possible. Advanced settings are hidden by default and only need to be opened when manual override is necessary.

                    **Typical outputs**
                    - Recovered DNA payload  
                    - Recovery summary  
                    - Decode check  
                    - Final restored file and preview, if successful
                    """
                )

        else:
            with st.container(border=True):
                st.subheader("📊 Recovery Analysis")
                st.markdown(
                    """
                    **Purpose**  
                    Evaluate whether recovery worked, how accurate it was, and whether the restored result can be trusted.

                    **How to use it**
                    1. Run recovery and decode in **NGS Decode**  
                    2. Open **Recovery Analysis**  
                    3. Review the summary, payload comparison, final decode result, and failure diagnosis  
                    4. Export reports if needed

                    **What this tab helps you answer**
                    - Was the DNA recovered successfully?  
                    - How similar is the recovered payload to the original payload?  
                    - Did zlib integrity pass?  
                    - Was a restored file produced?  
                    - Is the result trustworthy?

                    **Typical outputs**
                    - Recovery summary  
                    - Payload comparison  
                    - Final decode result  
                    - Failure diagnosis  
                    - Exportable analysis files
                    """
                )

    with tabs[2]:
        with st.container(border=True):
            st.subheader("Files needed in each tab")
            st.markdown(
                """
                **Encode**  
                - Original source file

                **NGS Prep**  
                - DNA payload from Encode  
                - or uploaded DNA text file

                **Sequencing / FASTQ**  
                - Fragment library from NGS Prep for simulation  
                - or external FASTQ files

                **NGS Decode**
                - For clean DNA decode: DNA payload text file  
                - For fragment-library decode: Manifest JSON from NGS Prep  
                - For FASTQ recovery: Manifest JSON from NGS Prep, R1 FASTQ, and optional R2 FASTQ
                """
            )

        with st.container(border=True):
            st.subheader("Where the files come from")
            st.markdown(
                """
                - **DNA payload text** comes from **Encode**  
                - **Manifest JSON**, **fragment FASTA**, and **fragment CSV** come from **NGS Prep**  
                - **R1/R2 FASTQ** and **FASTQ summary JSON** come from **Sequencing / FASTQ**  
                - **Recovered DNA payload**, **recovery report**, and **decode log** come from **NGS Decode** and **Recovery Analysis**
                """
            )

    with tabs[3]:
        with st.container(border=True):
            st.subheader("Common questions")
            st.markdown(
                """
                **Why does recovery fail even when FASTQ files exist?**  
                Recovery can fail if there are too many missing or corrupted fragments, if the protection mode is too weak for the error level, or if the reconstructed DNA does not pass the decode and zlib integrity checks.

                **Why do fragment statistics differ from the clean payload statistics?**  
                Fragment statistics are computed on full sequencing-ready fragments, not only on the clean payload. Primers, fragment identifiers, protection fields, and padding can change GC ratio and homopolymer behavior.

                **Why can one noisy run decode successfully while another fails?**  
                Recovery depends on the error profile, coverage, protection mode, and whether the reconstructed DNA remains valid enough to pass zlib integrity.

                **What is the best workflow for new users?**  
                Start with **Encode**, then **NGS Prep**, then **Sequencing / FASTQ**, then **NGS Decode**, and finally inspect the result in **Recovery Analysis**.
                """
            )

        with st.container(border=True):
            st.subheader("Recommended first test")
            st.markdown(
                """
                For a first test:
                1. Use a small file in **Encode**  
                2. Choose **None** or **CRC** in **NGS Prep**  
                3. Simulate FASTQ with light substitution noise in **Sequencing / FASTQ**  
                4. Recover and decode in **NGS Decode**  
                5. Inspect confidence and failure diagnosis in **Recovery Analysis**
                """
            )
