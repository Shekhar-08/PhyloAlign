import streamlit as st
from Bio import AlignIO, SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo import draw
from pathlib import Path
from io import BytesIO
from io import StringIO
import tempfile
import os
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import subprocess
import base64
from collections import Counter
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = 50_000_000

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="PhyloAlign", layout="wide", page_icon="🧬")

# Create a light green conceptual background image
fig, ax = plt.subplots(figsize=(10, 6), facecolor='#f0f9f0')  # Very light green background

# Add DNA sequence pattern (darker colors for contrast)
for i in range(15):
    y = np.linspace(0, 1, 100)
    x = np.sin(y * 20 + i) * 0.1 + i/15
    ax.plot(x, y, color='#1a5276', alpha=0.2, linewidth=2)  # Dark blue for contrast

# Add tree-like structure
for i in range(5):
    x = [0.5, 0.3 + i*0.1, 0.7 - i*0.1]
    y = [0.2, 0.5, 0.5]
    ax.plot(x, y, color='#27ae60', alpha=0.3, linewidth=3)  # Medium green

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
plt.tight_layout()
plt.savefig('bio_bg_lightgreen.png', transparent=True, dpi=300, facecolor=fig.get_facecolor())

def add_lightgreen_bg_image():
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("data:image/png;base64,{base64.b64encode(Path('bio_bg_lightgreen.png').read_bytes()).decode()}");
            background-size: cover;
            background-position: center;
            background-repeat: no-repeat;
            background-attachment: fixed;
            background-blend-mode: overlay;
            background-color: rgba(240, 249, 240, 0.95);
        }}
        .block-container {{
            background-color: rgba(255, 255, 255, 0.95);
            border-radius: 10px;
            padding: 2rem;
            border: 1px solid rgba(0, 0, 0, 0.1);
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.05);
        }}
        h1, h2, h3, h4, h5, h6 {{
            color: #1e8449 !important;
        }}
        p, .stText, .stMarkdown {{
            color: #333333 !important;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

# Call this function at the start of your app
add_lightgreen_bg_image()

# --- CLUSTALO PATH FOR LINUX ---
# Updated path for Linux environment
CLUSTALO_PATH = "clustalo"

# Verify Clustal Omega is executable
try:
    subprocess.run([CLUSTALO_PATH, "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except Exception as e:
    st.error(f"Clustal Omega not found or not executable at {CLUSTALO_PATH}. Please verify installation.")
    st.stop()

# --- SIDEBAR ---
st.sidebar.markdown("""
<style>
.sidebar-section {
    background-color: white;
    border-radius: 10px;
    padding: 10px;
    margin: 10px 0;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}
.sidebar-heading {
    color: black !important;
    font-weight: 800 !important;
    font-size: 18px !important;
    margin-bottom: 10px !important;
}
.sidebar-option {
    font-weight: 600 !important;
}
</style>
""", unsafe_allow_html=True)

# Main title
st.sidebar.markdown('<div class="sidebar-section"><h1 class="sidebar-heading">🧬 PhyloAlign v1.0</h1></div>', unsafe_allow_html=True)

# Menu section
with st.sidebar.expander("**NAVIGATION**", expanded=True):
    selection = st.radio(
        "Go to",
        [
            "Home",
            "Sequence Input",
            "Alignment Results",
            "Phylogenetic Tree",
            "About",
            "Feedback & Contact"
        ],
        key='menu_radio',
        label_visibility="collapsed"
    )

# Contacts section
st.sidebar.markdown("""
<div class="sidebar-section">
    <h3 class="sidebar-heading">CONTACTS</h3>
    <p class="sidebar-option">
        <a href="https://github.com/Shekhar-08" style="color: black; font-weight: 600;">GitHub</a> | 
        <a href="https://www.linkedin.com/in/shekhar-gudda-0299062ab/" style="color: black; font-weight: 600;">LinkedIn</a>
    </p>
</div>
""", unsafe_allow_html=True)

# Optional: Add a footer
st.sidebar.markdown("""
<div class="sidebar-section" style="text-align: center; padding: 5px;">
    <p style="font-weight: 600; margin: 0;">© 2025 PhyloAlign</p>
</div>
""", unsafe_allow_html=True)

# --- MAIN PAGE ---
st.markdown("""
    <div style='background-color:#E3F2FD; padding: 10px; border-radius: 10px; margin-bottom: 15px;'>
        <h1 style='color:#0D47A1; text-align: center;'>PhyloAlign</h1>
    </div>
""", unsafe_allow_html=True)

if selection == "Home":
    st.markdown("""
        <div style='background-color: #e6f4ea; padding: 20px; border-radius: 10px;'>
            <h3 style='color: #2c662d;'>Welcome to <strong>PhyloAlign</strong>: Multiple Sequence Alignment and Phylogenetic Analysis Web App</h3>
            <p>
                <strong>PhyloAlign</strong> is a Streamlit-based interactive web tool that enables researchers to perform Multiple Sequence Alignment (MSA) and construct Phylogenetic Trees using Clustal Omega and BioPython.
                This tool facilitates the alignment of DNA or protein sequences and visualizes both the aligned sequences and the evolutionary relationship among them.
            </p>
            <p>
                🎯 Key Features:
                <ul>
                    <li><strong>Manual or FASTA File for sequencesInput:</strong> 
                    You can input sequences either by pasting them manually into the text area or by uploading a FASTA file. Both options support standard FASTA formatting with headers and sequence lines.</li>
                    <li><strong>Clustal Omega Integration for alignment and guide tree generation:</strong> 
                    Clustal Omega is used to perform fast and accurate multiple sequence alignment. It also generates a guide tree to represent evolutionary relationships based on the input sequences.</li>
                    <li><strong>Color-coded alignment visualization:</strong> 
                    Aligned sequences are displayed using color codes to highlight nucleotide or amino acid similarities and differences. This helps in easily identifying conserved regions and mutations.</li>
                    <li><strong>Phylogenetic tree display:</strong> 
                    The phylogenetic tree visualizes the evolutionary relationships among the aligned sequences. It is generated using the guide tree output from Clustal Omega and displayed with BioPython.</li>
                </ul>
            </p>
        </div>
    """, unsafe_allow_html=True)

    st.markdown("---")
    
    with st.expander("📘 User Guide", expanded=False):
        st.markdown("""
        ### How to Use (3 Simple Steps)
        
        **1. Input Sequences**
        - **Option A: Paste Manually**
          1. Go to → *Sequence Input → Paste Sequences*
          2. Copy-paste in FASTA format:
             ```
             >Human_HBB
             ACGTACGT...
             ```
          3. Click *Submit Sequences*

        - **Option B: Upload File**
          1. Go to → *Sequence Input → Upload FASTA File*
          2. Select your `.fasta` file

        **2. View Alignment**
        - Automatic results in:
          - Color-coded sequences
          - Text alignment (FASTA)
          - Sequence logo (conservation)

        **3. Phylogenetic Tree**
        - Navigate to → *Phylogenetic Tree*
        - Interactive visualization appears
        - Download as Newick file

        ### Examples
        **Sample FASTA Format:**
        ```fasta
        >Human_Hemoglobin
        MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDL
        >Chimpanzee_Hemoglobin  
        MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDL
        ```

        ### Troubleshooting
        - **Format Error?** Ensure headers start with `>`
        - **Slow Performance?** Try <100 sequences
        - **Tree Not Loading?** Reduce sequence length
        - **Visualization Issues?** Download the full image
        """)

# --- SEQUENCE INPUT SECTION ---
if selection == "Sequence Input":
    st.header("🧬 Enter or Upload Sequences")
   # Example FASTA sequences
    example_sequences = """
>Human_HBB
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
>Chimpanzee_HBB
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
>Gorilla_HBB
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
>Mouse_HBB
MVHLTDAEKAAVNGLWGKVNPDDVGGEALGRLLVVYPWTQRYFDSFGDLSSASAVMGNPK""" 
   
    input_method = st.radio("Choose input method:", ["Paste Sequences", "Upload FASTA File"])

    sequences = []
    if input_method == "Paste Sequences":
        with st.expander("💡 Example FASTA Format (Click to view/copy)"):
            st.code(example_sequences, language="text")
            
        seq_input = st.text_area("Paste sequences in FASTA format:", 
                                height=200,
                                placeholder="Paste your sequences here in FASTA format...\n\nExample:\n>Sequence1\nACGTACGTACGT\n>Sequence2\nTGCAACGTACGT")
        if st.button("Submit Sequences") and seq_input:
            try:
                sequences = list(SeqIO.parse(StringIO(seq_input), "fasta"))
                st.success(f"✅ Loaded {len(sequences)} sequences.")
                st.session_state["sequences"] = sequences
                # Clear previous results when new sequences are loaded
                if "alignment" in st.session_state:
                    del st.session_state["alignment"]
                if "tree" in st.session_state:
                    del st.session_state["tree"]
            except Exception as e:
                st.error(f"Error parsing sequences: {e}")

    elif input_method == "Upload FASTA File":
        uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa", "txt"])
        if uploaded_file is not None:
            sequences = list(SeqIO.parse(uploaded_file, "fasta"))
            st.success(f"✅ Loaded {len(sequences)} sequences.")
            st.session_state["sequences"] = sequences
            # Clear previous results when new sequences are loaded
            if "alignment" in st.session_state:
                del st.session_state["alignment"]
            if "tree" in st.session_state:
                del st.session_state["tree"]

    if sequences:
        st.subheader("Input Sequences Preview:")
        st.code("".join(f">{record.id}\n{record.seq}\n" for record in sequences[:5]) + 
               ("\n...\n" if len(sequences) > 5 else ""), language="text")

# Function to run Clustal Omega and store results
def run_clustal_omega(sequences):
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = os.path.join(tmpdir, "input.fasta")
        output_path = os.path.join(tmpdir, "aligned.fasta")
        tree_path = os.path.join(tmpdir, "tree.newick")

        with open(input_path, "w") as f:
            SeqIO.write(sequences, f, "fasta")

        # Use absolute paths for Linux compatibility
        cmd = [
            CLUSTALO_PATH,
            "-i", os.path.abspath(input_path),
            "-o", os.path.abspath(output_path),
            "--guidetree-out", os.path.abspath(tree_path),
            "--force",
            "--auto",
            "--verbose"
        ]

        try:
            # For Linux, we need to ensure proper permissions
            os.chmod(input_path, 0o644)
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                st.error(f"Clustal Omega error:\n{result.stderr}")
                return None, None
            
            alignment = AlignIO.read(output_path, "fasta")
            tree = Phylo.read(tree_path, "newick")
            
            return alignment, tree
            
        except Exception as e:
            st.error(f"❌ Error running Clustal Omega: {e}")
            return None, None

# --- ALIGNMENT RESULTS ---
if selection == "Alignment Results":
    if "sequences" not in st.session_state:
        st.warning("⚠️ Please input sequences first on the Sequence Input page.")
        st.stop()
    
    sequences = st.session_state["sequences"]
    
    if len(sequences) < 2:
        st.error("❌ You need at least 2 sequences for alignment.")
        st.stop()
    
    if "alignment" not in st.session_state or "tree" not in st.session_state:
        with st.spinner("Running Clustal Omega for alignment and tree generation..."):
            alignment, tree = run_clustal_omega(sequences)
            if alignment is not None and tree is not None:
                st.session_state["alignment"] = alignment
                st.session_state["tree"] = tree
                st.success("✅ Alignment and tree generated successfully!")
            else:
                st.error("Failed to generate alignment and tree.")
                st.stop()
    
    alignment = st.session_state["alignment"]
    
    st.subheader("🧬 Aligned Sequences:")
    output_buffer = StringIO()
    AlignIO.write(alignment, output_buffer, "fasta")
    st.code(output_buffer.getvalue(), language="text")
    
    # Download button for aligned sequences
    st.download_button(
        label="📥 Download Aligned FASTA",
        data=output_buffer.getvalue(),
        file_name="aligned.fasta",
        mime="text/plain"
    )
    
    # Color-coded alignment visualization
    st.subheader("🎨 Color-coded Alignment")
    
    def plot_alignment(alignment):
        color_scheme = {
            'A': '#1E90FF', 'C': '#FF8C00', 'G': '#32CD32', 'T': '#DC143C', 'U': '#DC143C',
            'N': '#A9A9A9', '-': '#FFFFFF'
        }
        
        MAX_SEQS = 100
        MAX_BASES = 500
        
        n_seqs = min(len(alignment), MAX_SEQS)
        seq_len = min(len(alignment[0]), MAX_BASES)
        
        fig, ax = plt.subplots(figsize=(min(seq_len*0.2, 100), min(n_seqs*0.5, 50)), dpi=100)
        ax.set_xlim(0, seq_len)
        ax.set_ylim(0, n_seqs)
        ax.axis('off')
        
        for y, record in enumerate(alignment[:MAX_SEQS]):
            for x, char in enumerate(str(record.seq[:MAX_BASES])):
                color = color_scheme.get(char.upper(), '#000000')
                ax.add_patch(plt.Rectangle((x, n_seqs-y-1), 1, 1, color=color, linewidth=0))
                if n_seqs <= 30 and seq_len <= 100:
                    ax.text(x+0.5, n_seqs-y-0.5, char, ha='center', va='center',
                           fontsize=8, fontfamily='monospace',
                           color='white' if char.upper() in ['A','C','G','T','U'] else 'black')
        
        if len(alignment) > MAX_SEQS or len(alignment[0]) > MAX_BASES:
            ax.text(0.5, 0.5, 
                   f"Showing first {MAX_SEQS} of {len(alignment)} seqs\n"
                   f"First {MAX_BASES} of {len(alignment[0])} bases",
                   ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        return fig

    # Display the visualization
    if "alignment" in st.session_state:
        alignment = st.session_state["alignment"]
        
        if len(alignment) > 100 or len(alignment[0]) > 500:
            st.warning("Large alignment detected - showing simplified view")
        
        try:
            with st.spinner("Generating visualization..."):
                fig = plot_alignment(alignment)
                st.pyplot(fig, use_container_width=True)
                
                # Add download button for color-coded alignment
                buf = BytesIO()
                fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
                plt.close(fig)
                st.download_button(
                    label="📥 Download Color-coded Alignment (PNG)",
                    data=buf.getvalue(),
                    file_name="color_coded_alignment.png",
                    mime="image/png"
                )
        except Exception as e:
            st.error(f"Visualization error: {str(e)}")
            st.info("Showing sequence logo instead...")
            try:
                counts_df = logomaker.alignment_to_matrix([str(rec.seq[:200]) for rec in alignment[:50]])
                logo_fig, ax = plt.subplots(figsize=(20, 5))
                logomaker.Logo(counts_df, ax=ax)
                st.pyplot(logo_fig)
                
                # Add download button for sequence logo
                buf = BytesIO()
                logo_fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
                plt.close(logo_fig)
                st.download_button(
                    label="📥 Download Sequence Logo (PNG)",
                    data=buf.getvalue(),
                    file_name="sequence_logo.png",
                    mime="image/png"
                )
            except Exception as e2:
                st.error(f"Could not generate any visualization: {str(e2)}")
        
        st.markdown("---")
    
    # Sequence Logo
    def create_scrollable_logo(alignment, max_positions=500, logo_height=8):
        """Generate a scrollable sequence logo with increased height"""
        # Calculate dimensions
        seq_length = min(len(alignment[0]), max_positions)
        logo_width = max(50, seq_length * 0.3)  # Wider than before
    
        # Create figure with increased height
        fig, ax = plt.subplots(figsize=(logo_width, logo_height), dpi=100)
    
        try:
            # Create frequency matrix
            counts_df = logomaker.alignment_to_matrix(
                [str(rec.seq[:seq_length]) for rec in alignment]
            )
        
            # Generate taller logo
            logo = logomaker.Logo(counts_df,
                                ax=ax,
                                font_name='Arial',
                                color_scheme='chemistry',
                                vpad=0.2,  # More vertical padding
                                width=0.9,
                                flip_below=False)
        
            # Enhanced styling
            logo.style_spines(visible=False)
            logo.ax.set_ylabel('Bits', fontsize=14, fontweight='bold')
            logo.ax.set_xlabel('Position', fontsize=12, labelpad=10)
            plt.title('Sequence Conservation', fontsize=16, pad=20)
        
            # Adjust layout
            plt.tight_layout()
            return fig
        
        except Exception as e:
            plt.close()
            raise e

    # In your Streamlit app:
    st.subheader("📊 Interactive Sequence Logo")

    if "alignment" in st.session_state:
        alignment = st.session_state["alignment"]
    
        # Warning for large alignments
        if len(alignment[0]) > 500:
            st.warning(f"Displaying first 500 of {len(alignment[0])} positions. Download for full analysis.")
    
        try:
            with st.spinner("Generating enhanced visualization..."):
                fig = create_scrollable_logo(alignment, logo_height=10)  # Increased height
            
                # Create scrollable container with custom CSS
                st.markdown("""
                <style>
                .logo-container {
                    overflow-x: auto;
                    border: 1px solid #e1e4e8;
                    border-radius: 8px;
                    padding: 15px;
                    background: white;
                    margin: 10px 0;
                }
                .logo-container img {
                    min-width: 1500px;  # Ensures scrolling
                    height: auto;
                    display: block;
                }
                </style>
                """, unsafe_allow_html=True)
            
                with st.container():
                    st.markdown('<div class="logo-container">', unsafe_allow_html=True)
                    st.pyplot(fig, use_container_width=False)
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                # Add download button for sequence logo
                buf = BytesIO()
                fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
                plt.close(fig)
                st.download_button(
                    label="📥 Download Sequence Logo (PNG)",
                    data=buf.getvalue(),
                    file_name="sequence_logo.png",
                    mime="image/png"
                )
                
        except Exception as e:
            st.error(f"Visualization error: {str(e)}")
        
        finally:
            plt.close('all')

# --- PHYLOGENETIC TREE ---
if selection == "Phylogenetic Tree":
    if "tree" not in st.session_state:
        st.warning("⚠️ Please generate an alignment first on the Alignment Results page.")
        st.stop()
    
    tree = st.session_state["tree"]
    
    st.subheader("🌳 Phylogenetic Tree")
    
    # Customize tree appearance
    tree.ladderize()  # Put more closely related sequences on the same branch
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_title("Phylogenetic Tree", pad=20)
    
    # Draw the tree
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Customize appearance
    ax.set_facecolor('#f5f5f5')
    fig.patch.set_facecolor('#f5f5f5')
    
    st.pyplot(fig)
    
    # Download button for tree
    tree_buffer = StringIO()
    Phylo.write(tree, tree_buffer, "newick")
    st.download_button(
        label="📥 Download Tree (Newick format)",
        data=tree_buffer.getvalue(),
        file_name="phylogenetic_tree.newick",
        mime="text/plain"
    )

# --- ABOUT ---
if selection == "About":
    st.subheader("About PhyloAlign")
    st.markdown("""
        PhyloAlign is an intuitive bioinformatics web tool designed to simplify multiple sequence alignment (MSA) and phylogenetic analysis. Built with Streamlit and Biopython, it leverages Clustal Omega to enable researchers, educators, and students to align DNA/protein sequences and visualize evolutionary relationships through interactive trees and color-coded alignments—all without command-line expertise. The tool offers instant visual feedback via sequence logos and publication-ready graphics, making it ideal for studying mutations, conserved domains, or teaching core bioinformatics concepts. By combining accuracy with accessibility, PhyloAlign bridges the gap between advanced computational biology and user-friendly exploration.
        
        **Version:** 1.0  
        **Author:** Shekhar Gudda  
        
        Built with:
        - Streamlit
        - Biopython
        - Clustal Omega
        - Matplotlib
    """)
    
    st.markdown("""
        <h3 style='color:#27ae60;'>🔬 Tool Overview</h3>
        <p>
            PhyloAlign is a bioinformatics web application designed for biological sequence analysis, offering:
        </p>
        <ul>
            <li><strong>Multiple Sequence Alignment (MSA)</strong> using Clustal Omega algorithm</li>
            <li><strong>Phylogenetic Tree Visualization</strong> from alignment results</li>
            <li><strong>Interactive Visualizations</strong> including color-coded alignments and sequence logos</li>
            <li><strong>User-Friendly Interface</strong> for researchers at all bioinformatics skill levels</li>
        </ul>
        </p>
        <h3 style='color:#27ae60;'>🧬 Scientific Applications</h3>
        <p>
            Ideal for:
        </p>
        <ul>
            <li>Evolutionary biology studies</li>
            <li>Conserved domain identification</li>
            <li>Primer design validation</li>
            <li>Protein/DNA sequence comparison</li>
            <li>Teaching molecular biology concepts</li>
        </ul>
        </p>
        <h3 style='color:#27ae60;'>⚙️ Technical Specifications</h3>
        <ul>
            <li><strong>Backend:</strong> Python 3.9+ with Biopython</li>
            <li><strong>Alignment Engine:</strong> Clustal Omega (v1.2.2)</li>
            <li><strong>Visualization:</strong> Matplotlib, Logomaker</li>
            <li><strong>Web Framework:</strong> Streamlit</li>
        </ul>
        </p>
        <div style='background-color:#e8f5e9; padding: 15px; border-radius: 8px; margin-top: 20px;'>
            <h3 style='color:#1e8449;'>👨‍💻 About the Author</h3>
            <div style='display: flex; align-items: center; gap: 20px; margin-bottom: 15px;'>
                <img src='https://avatars.githubusercontent.com/u/210685307?s=400&u=116167314105ea4825f635fe91de751cc2cda4f5&v=4' alt='Shekhar Gudda' style='width: 150px; height: 150px; border-radius: 50%; object-fit: cover; border: 3px solid #27ae60;'>
                <div>
                    <p style='margin: 0; font-size: 1.2em; font-weight: bold;'>Shekhar Gudda</p>
                    <p style='margin: 0; color: #555;'>Bioinformatics Researcher | Computational Biologist</p>
                </div>
            </div>
                <p>I am deeply passionate about creating accessible, intelligent tools to bridge the gap between computational analysis and biological discovery. PhyloAlign is a web-based application developed as part of my ongoing learning and contribution to the bioinformatics community. This tool enables users to perform multiple sequence alignment (MSA), visualize conserved regions, and generate phylogenetic trees—all through a user-friendly interface powered by Streamlit and Biopython.</p>
                <p>The motivation behind building PhyloAlign stems from the challenges I observed during genomic data interpretation and evolutionary analysis. My goal is to help students, educators, and researchers efficiently explore sequence relationships and evolutionary patterns.</p>
                </p>
                <strong>Education:</strong> Curremtly pursing Master's degree in Bioinformatics<br>
                <strong>Institution:</strong> DES Pune University<br>
                <strong>Research Focus:</strong> Evolutionary genomics, Sequence analysis, Phylogenetics<br><br>
                </p>
                <strong>Contact:</strong> 
                <a href="mailto:shekhargudda844@gmail.com">shekhargudda844@gmail.com</a><br>
                <strong>GitHub:</strong> <a href="https://github.com/Shekhar-08">github.com/Shekhar-08</a><br>
                <strong>LinkedIn:</strong> <a href="https://www.linkedin.com/in/shekhar-gudda-0299062ab/">https://www.linkedin.com/in/shekhar-gudda-0299062ab/</a>
            </p>
        </div>
        </p>
                <hr style='border: 0.5px solid #27ae60; margin: 20px 0; opacity: 0.3;'>
                </p>
                <h4 style='color:#1e8449; margin-top: 0;'>🙏 Acknowledgment</h4>
                <p>
                    I would like to sincerely thank <strong>Dr. Kushagra Kashyap</strong>, Assistant Professor, 
                    School of Science and Mathematics, DES Pune University, for his valuable guidance, 
                    support, and encouragement throughout the development of PhyloAlign. His insights in 
                    bioinformatics and computational biology have been instrumental in shaping this tool.
                    <br><br>
                    <strong>LinkedIn:</strong> <a href="https://linkedin.com/in/dr-kushagra-kashyap-b230a3bb" target="_blank">linkedin.com/in/dr-kushagra-kashyap</a>
                </p>
                <p>
                    Additional thanks to:
                </p>
                <ul style='margin-top: 5px;'>
                    <li>The open-source bioinformatics community for their invaluable tools</li>
                    <li>Streamlit and Biopython developers for making this technology accessible</li>
                </ul>
                <p style='font-size: 0.9em; margin-bottom: 0;'>
                    <em>This project was developed as part of my Master's research in Bioinformatics at DES Pune University.</em>
                </p>
            </div>
            <hr style='border: 0.5px solid #27ae60; margin: 20px 0; opacity: 0.3;'>
            </p>
        <div style='margin-top: 20px; font-size: 0.9em; color: #555;'>
            <p>
                <strong>Version:</strong> 1.0 | <strong>Last Updated:</strong> May 2025<br>
            </p>
        </div>
    </div>
    """, unsafe_allow_html=True)

# --- FEEDBACK SECTION ---
if selection == "Feedback & Contact":
    st.subheader("📬 Feedback & Contact")
    st.markdown("""
        We welcome your feedback to improve PhyloAlign!
        
        Please report any issues or suggestions to:
        - Email: shekhargudda844@gmail.com
        - GitHub: (https://github.com/Shekhar-08)
        - LinkedIn: (https://www.linkedin.com/in/shekhar-gudda-0299062ab/)
    """)
    
    with st.form("feedback_form"):
        name = st.text_input("Your Name")
        email = st.text_input("Your Email")
        feedback = st.text_area("Your Feedback or Questions")
        submitted = st.form_submit_button("Submit Feedback")
        if submitted:
            st.success("Thank you for your feedback! We'll get back to you soon.")
            

