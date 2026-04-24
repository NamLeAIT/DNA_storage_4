import streamlit as st

def render_software():
    st.title("Software & Resources")
    st.caption("Offline-first tooling for cleanroom and internet-restricted environments.")

    dl = st.container(border=True)
    with dl:
        st.subheader("💾 Download Software")
        st.write(
            """
            Our software is developed in-house and provided for offline use, enabling reliable operation
            in internet-restricted environments such as cleanrooms, or whenever users require stable,
            uninterrupted performance independent of network connectivity.
            """
        )
        col_v1, col_v2, col_v3 = st.columns(3, gap="large")

        with col_v1:
            st.markdown("##### DDSS ver.01")
            st.caption("Release Date: 2025.11")
            st.button("📥 Download v01", key="dl_v1", help="Stable version for basic encoding", use_container_width=True)

        with col_v2:
            st.markdown("##### DDSS ver.02")
            st.caption("Release Date: 2025.12")
            st.button("📥 Download v02", key="dl_v2", help="Added Reed-Solomon support", use_container_width=True)

        with col_v3:
            st.markdown("##### DDSS ver.03")
            st.caption("Latest Version")
            st.button("🚀 Download v03", key="dl_v3", help="Optimized for large video files", use_container_width=True)

        st.write("If you need assistance or encounter issues, please contact the lab team.")

    st.divider()

    tools = st.container(border=True)
    with tools:
        st.subheader("🔗 Related Tools")
        st.caption("Supporting tools and platforms for development, reproducibility, and bioinformatics.")
        link_col1, link_col2, link_col3 = st.columns(3, gap="large")

        with link_col1:
            st.markdown(
                """
                **Programming & Environment**
                - [Python Official](https://www.python.org/)
                - [Anaconda Distribution](https://www.anaconda.com/)
                - [Visual Studio Code](https://code.visualstudio.com/)
                """
            )

        with link_col2:
            st.markdown(
                """
                **AI & Support**
                - [ChatGPT (OpenAI)](https://chatgpt.com/)
                - [Claude AI](https://claude.ai/)
                - [GitHub Copilot](https://github.com/features/copilot)
                """
            )

        with link_col3:
            st.markdown(
                """
                **Bioinformatics Resources**
                - [Biopython Project](https://biopython.org/)
                - [NCBI Database](https://www.ncbi.nlm.nih.gov/)
                - [Hugging Face](https://huggingface.co/)
                """
            )

        st.button("Check for Updates", use_container_width=True)
