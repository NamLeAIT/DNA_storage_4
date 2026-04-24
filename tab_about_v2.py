
import streamlit as st


def render_about():
    st.title("About the DNA Data Storage Platform")
    st.caption("A sequencing-aware DNA data storage workflow for encoding, fragment preparation, sequencing simulation, recovery, and post-recovery analysis.")

    with st.container(border=True):
        st.subheader("Platform Overview")
        st.markdown(
            """
            This platform is designed as a practical DNA data storage workflow simulator that converts digital files
            into DNA payloads, prepares sequencing-ready fragment libraries, simulates sequencing-style FASTQ reads,
            reconstructs DNA from noisy inputs, and evaluates recovery reliability.

            The current system emphasizes the full computational chain from:
            **source file → encoded DNA payload → NGS-ready fragment library → FASTQ reads → recovered DNA → restored file**.

            In addition to baseline sequencing tests, the platform also supports protection-based recovery modes such as
            **CRC**, **Parity Fragments**, and **Reed–Solomon**, allowing users to study both recovery robustness and the
            storage overhead introduced by stronger protection schemes.
            """
        )

    with st.container(border=True):
        st.subheader("What makes this version different")
        st.markdown(
            """
            Compared with earlier versions, this version provides:

            - A clearer **pipeline-based interface**  
            - A dedicated **Sequencing / FASTQ** tab for sequencing simulation and FASTQ import  
            - A dedicated **Recovery Analysis** tab for post-recovery evaluation  
            - Natural-language labels and cleaner UI wording  
            - Hidden advanced settings by default, with automatic reuse of stored parameters where appropriate  
            - Preview support for recovered images, audio, video, and text files
            """
        )

    with st.container(border=True):
        st.subheader("Main workflow supported by the platform")
        st.markdown(
            """
            1. **Encode** — compress and convert a source file into a DNA payload  
            2. **NGS Prep** — split the DNA payload into sequencing-ready fragments and apply a protection mode  
            3. **Sequencing / FASTQ** — simulate noisy sequencing reads or import external FASTQ files  
            4. **NGS Decode** — recover DNA and optionally restore the original file  
            5. **Recovery Analysis** — inspect the quality, confidence, and likely failure points of the recovered result
            """
        )

    with st.container(border=True):
        st.subheader("Members and Contact")
        st.write("If you have questions, suggestions, or feedback, please contact:")

        c1, c2 = st.columns(2, gap="large")
        with c1:
            st.markdown(
                """
                **Dr. Nguyen Kim Uyen**  
                Research Developer  
                Email: kimuyendlu@gmail.com
                """
            )
        with c2:
            st.markdown(
                """
                **Mr Le Quoc Nam**  
                Research Support  
                Email: —
                """
            )

    with st.container(border=True):
        st.subheader("Affiliation and Location")
        st.markdown(
            """
            **Affiliation:** Department of Physics, Sungkyunkwan University  
            **Address:** 2066, Seobu-ro, Jangan-gu, Suwon, Gyeonggi-do, 16419, Republic of Korea
            """
        )

        map_html = """
<iframe src="https://www.google.com/maps/embed?pb=!1m18!1m12!1m3!1d3171.303975005891!2d126.97210167637823!3d37.2938883394593!2m3!1f0!2f0!3f0!3m2!1i1024!2i768!4f13.1!3m3!1m2!1s0x357b56b21761867f%3A0xb38ea754e92d9bb0!2sSungkyunkwan%20University%20(Natural%20Sciences%20Campus)!5e0!3m2!1sen!2skr!4v1700000000000!5m2!1sen!2skr"
    width="100%"
    height="430"
    style="border:0; border-radius: 14px;"
    allowfullscreen=""
    loading="lazy"
    referrerpolicy="no-referrer-when-downgrade">
</iframe>
"""
        st.components.v1.html(map_html, height=480)
