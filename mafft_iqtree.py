import streamlit as st
import subprocess
import tempfile
import os
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Align.Applications import MafftCommandline


# File paths for MAFFT and IQTREE
mafft_path = "GUIMAFFT/mafft.bat"
iqtree_path = "GUIMAFFT/iqtree"

# Formatting
st.markdown(
    """
    <link href='https://fonts.googleapis.com/css?family=Oxanium' rel='stylesheet'>

    <style>
    .stApp {
        background-size: cover;
        background-position: center;
        background-repeat: no-repeat;
        background-attachment: fixed;
    }
    .title-panel {
        background: linear-gradient(90deg, rgba(160,237,255,1) 15%, rgba(115,222,81,1) 62%, rgba(0,212,255,1) 100%);
        padding: 10px;
        border-radius: 8px;
        text-align: center;
        font-size: 48px;
        color: #fff;
        margin-bottom: 5px;
        font-family: 'Oxanium';
    }
    </style>
    <div class="title-panel">
        <h1>Alignment and Phylogeny</h1> <p>MAFFT and IQ-TREE</p>
    </div>
    """,
    unsafe_allow_html=True
)

uploaded_files = st.file_uploader("Upload your FASTA file(s)", type=["fasta", "fa"], accept_multiple_files=True)

algorithm_options = {
    "--auto": {},
    "FFT-NS-1 (fast)": {"retree": 1},
    "FFT-NS-2 (default)": {"retree": 2},
    "G-INS-i (accurate)": {"globalpair": True, "maxiterate": 16},
    "L-INS-i (accurate)": {"localpair": True, "maxiterate": 16},
    "E-INS-i (accurate)": {"genafpair": True, "maxiterate": 16}
}

if uploaded_files:
    choice = st.radio("Choose alignment algorithm and click 'Run MAFFT Alignment'", list(algorithm_options.keys()),
                      help="Choose the algorithm for MAFFT.")
    selected_options = algorithm_options.get(choice, None)
    run_button = st.button("Run MAFFT Alignment")

    if run_button and selected_options is not None:
        for uploaded_file in uploaded_files:
            # Save the uploaded file to a temporary location
            with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_fasta:
                temp_fasta.write(uploaded_file.read())
                fasta_path = temp_fasta.name

            output_path = tempfile.mktemp(suffix=".aligned.fasta")
            status_message = st.empty()
            status_message.info(f"Running MAFFT alignment on {uploaded_file.name}...")

            # Set up MAFFT command line with selected options
            mafft_cline = MafftCommandline(input=fasta_path, **selected_options)

            try:
                # Execute MAFFT and capture output
                stdout, stderr = mafft_cline()
                with open(output_path, 'w') as aligned_file:
                    aligned_file.write(stdout)

                status_message.success(f"MAFFT alignment completed for {uploaded_file.name}!", icon="✅")

                # Store alignment in session state to prevent reset
                st.session_state[f"aligned_{uploaded_file.name}"] = stdout

                # Display alignment results and download button
                st.sidebar.subheader(f"Alignment Results for {uploaded_file.name}")
                st.sidebar.text(st.session_state[f"aligned_{uploaded_file.name}"])

                st.download_button(f"Download Aligned {uploaded_file.name}",
                                   data=st.session_state[f"aligned_{uploaded_file.name}"],
                                   file_name=f"aligned_{uploaded_file.name}", mime="text/fasta")
            except Exception as e:
                status_message.error(f"MAFFT alignment failed for {uploaded_file.name}: {str(e)}")


# IQ-TREE options
if uploaded_files and any(f"aligned_{file.name}" in st.session_state for file in uploaded_files):
    iqtree_options = st.radio("Choose IQ-TREE model", ["Auto", "GTR+G", "HKY", "JC"],
                              help="Choose the model for phylogenetic tree construction.")

    run_iqtree_button = st.button("Construct Phylogenetic Trees with IQ-TREE")

    if run_iqtree_button:
        for uploaded_file in uploaded_files:
            alignment_content = st.session_state.get(f"aligned_{uploaded_file.name}", None)

            if alignment_content:
                # Status message for IQ-TREE running
                status_message = st.empty()
                status_message.info(f"Running IQ-TREE for {uploaded_file.name}...")

                # Save alignment data to a temporary file
                with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_fasta:
                    temp_fasta.write(alignment_content.encode())
                    fasta_path_iqtree = temp_fasta.name

                output_tree_path = tempfile.mktemp(suffix=".treefile")

                # IQ-TREE command
                model_option = "-m MFP" if iqtree_options == "Auto" else f"-m {iqtree_options}"
                iqtree_command = f"{iqtree_path} -s {fasta_path_iqtree} {model_option} -pre {output_tree_path}"

                # Run IQ-TREE and capture output
                process = subprocess.Popen(iqtree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()

                if process.returncode != 0:
                    st.error(f"IQ-TREE failed with error: {stderr.decode()}")
                else:
                    status_message.success("IQ-TREE tree construction completed!", icon="✅")

                    try:
                        # Check if output tree file exists
                        if os.path.exists(f"{output_tree_path}.treefile"):
                            tree = Phylo.read(f"{output_tree_path}.treefile", "newick")
                            tree.ladderize()
                            fig, ax = plt.subplots(figsize=(10, 8))


                            Phylo.draw(tree, do_show=False, axes=ax)



                            st.subheader(f"Phylogenetic Tree for {uploaded_file.name}")
                            st.pyplot(fig)

                            # Provide download options for the tree and PHYLIP files
                            st.download_button(f"Download Tree File ({uploaded_file.name})",
                                               data=open(f"{output_tree_path}.treefile").read(),
                                               file_name=f"{uploaded_file.name}.treefile", mime="text/plain")

                            phylip_file = f"{output_tree_path}.phy"
                            if os.path.exists(phylip_file):
                                st.download_button(f"Download PHYLIP File ({uploaded_file.name})",
                                                   data=open(phylip_file).read(),
                                                   file_name=f"{uploaded_file.name}.phy", mime="text/plain")
                        else:
                            st.error("Tree file not found after IQ-TREE execution.")
                    except Exception as e:
                        st.error(f"Failed to load the tree: {e}")

                # Clean up temporary files
                os.remove(fasta_path_iqtree)
                if os.path.exists(output_tree_path):
                    os.remove(output_tree_path)
                if os.path.exists(f"{output_tree_path}.treefile"):
                    os.remove(f"{output_tree_path}.treefile")
                if os.path.exists(f"{output_tree_path}.phy"):
                    os.remove(f"{output_tree_path}.phy")
