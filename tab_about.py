import streamlit as st

def render_about():
    st.title("About the DNA Data Storage Framework")
    st.caption("An end-to-end DNA data storage framework integrating compression, headerless zlib framing, DNA sequence design, and recovery validation.")

    overview = st.container(border=True)
    with overview:
        st.subheader("System Overview")
        st.markdown(
            """
            This framework is designed as an end-to-end DNA data storage system that transforms digital files into DNA sequences
            through data-type-aware compression, headerless zlib framing, and biologically constrained binary-to-DNA encoding.

            The system emphasizes the full computational chain from raw input to recoverable DNA payload, including automatic
            benchmark-based representation selection, payload construction without custom in-band headers, sequence generation
            under GC and homopolymer considerations, and decoding with integrity-aware restoration.
            """
        )

    st.divider()

    members = st.container(border=True)
    with members:
        st.subheader("Members and Contact")
        st.write("If you have questions or feedback, please contact:")

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

    st.divider()

    loc = st.container(border=True)
    with loc:
        st.subheader("Location")
        st.markdown(
            """
            **Address:** Department of Physics, Sungkyunkwan University  
            2066, Seobu-ro, Jangan-gu, Suwon, Gyeonggi-do, 16419, Republic of Korea
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
