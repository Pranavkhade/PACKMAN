#!/usr/bin/env python
"""The PACKMAN Streamlit Web Application.

This is a modern web interface for PACKMAN functionality using Streamlit.

Features:
    * Hinge Prediction with ScanAlpha
    * hdANM molecular dynamics
    * Voronoi Packing Entropy calculation
    * Interactive 3D protein visualization

Authors:
    * Pranav Khade (https://github.com/Pranavkhade)
"""

import streamlit as st
import numpy as np
import os
import tempfile
import logging
from pathlib import Path
from io import StringIO

# PACKMAN imports
from packman import molecule
from packman.apps import predict_hinge
from packman.anm import hdANM
from packman.entropy import PackingEntropy

# For visualization
import py3Dmol
import stmol

# Set page config
st.set_page_config(
    page_title="PACKMAN Web Interface",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Global variables
hinge_colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple', 'cyan', 'magenta']

# Get PACKMAN version
def get_version():
    """Get PACKMAN version"""
    try:
        import packman
        return packman.__version__
    except:
        return "Unknown"

# Custom CSS
st.markdown("""
    <style>
    .main {
        padding: 0rem 0rem;
    }
    .stTabs [data-baseweb="tab-list"] button [data-testid="stMarkdownContainer"] p {
        font-size: 1.25rem;
    }
    </style>
    """, unsafe_allow_html=True)

def render_protein_with_hinges(pdb_path, hinge_residues, chain_id=None):
    """Render protein structure with highlighted hinges"""
    
    if( os.path.splitext(pdb_path)[-1] == '' ):
        pdb_path = pdb_path + '.cif'
    else:
        pdb_path = pdb_path
        

    extension = os.path.splitext(pdb_path)[1].lower()

    if extension not in ['.pdb', '.cif']:
        st.error("Unsupported file format for visualization. Please use PDB or CIF files.")
        return
    
    try:
        if extension == '.pdb':
            view = py3Dmol.view(width=800, height=600)
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
            view.addModel(pdb_content, 'pdb')
        else:
            view = py3Dmol.view(width=800, height=600)
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
            view.addModel(pdb_content, 'cif')
        
        # Color protein white
        view.setStyle({'ss': 'h'}, {'cartoon': {'color': 'white'}})
        view.setStyle({'ss': 's'}, {'cartoon': {'color': 'white'}})
        view.setStyle({'ss': 'c'}, {'cartoon': {'color': 'white'}})
        
        # Color hinges
        for idx, (start, end) in enumerate(hinge_residues):
            color = hinge_colors[idx % len(hinge_colors)]
            for res_id in range(start, end + 1):
                if chain_id:
                    view.setStyle({'chain': chain_id, 'resi': res_id}, {'cartoon': {'color': color}})
                else:
                    view.setStyle({'resi': res_id}, {'cartoon': {'color': color}})
        
        view.zoomTo()
        stmol.showmol(view, height=600, width=800)
        
    except Exception as e:
        st.error(f"Error rendering structure: {str(e)}")

def render_generated_movie(cif_file):
    """Render generated movie from hdANM"""
    
    mol = molecule.load_cif(cif_file)

    # save individual frames as cif files and add slider to select frame
    all_files = []
    for numi,i in enumerate(*mol.get_models()):
        frame = molecule.Protein(str(numi), [i])
        frame.write_pdb(f"frame_{numi}.pdb")
        all_files.append(f"frame_{numi}.pdb")

    try:
        view = py3Dmol.view(width=800, height=600)

        all_models = ""
        for numf, f in enumerate(all_files):
            with open(f, 'r') as file:
                local_file_content = file.read()
                local_file_content = local_file_content.replace("Model\t0\n", f"MODEL {numf}\n")
                all_models += local_file_content
                #st.info(local_file_content.split('\n')[0])
        all_models += "END\n"  # End of PDB file

        # save file
        with open("movie_frames.pdb", 'w') as f:
            f.write(all_models)
    
        view.addModelsAsFrames(all_models, 'pdb')
        
        view.setStyle({}, {'cartoon': {'color': 'spectrum'}})

        view.zoomTo()
        view.animate({'loop': 'forward', 'interval': 50})
        stmol.showmol(view, height=600, width=800)

    except Exception as e:
        st.error(f"Error rendering movie: {str(e)}")


def scan_alpha_hinge_prediction(backbone_atoms, alpha_start, alpha_stop, step_size):
    """
    Scan alpha values for hinge prediction.
    Uses the ScanAlpha algorithm to identify robust hinges.
    """
    chain_id = backbone_atoms[0].get_parent().get_parent().get_id()
    
    # Scan through alpha values
    for alpha_val in np.arange(alpha_start, alpha_stop, step_size):
        alpha_val = np.around(alpha_val, decimals=1)
        try:
            predict_hinge(backbone_atoms, Alpha=alpha_val, outputfile=open(os.devnull, 'w'))
        except:
            continue
    
    hinges = backbone_atoms[0].get_parent().get_parent().get_hinges()
    
    # Remove insignificant hinges with p-value < 0.05
    pre_overlap = []
    for i in range(0, len(hinges)):
        for j in range(i + 1, len(hinges)):
            h1 = [k.get_id() for k in hinges[i].get_elements()]
            h2 = [k.get_id() for k in hinges[j].get_elements()]
            if (hinges[i].get_pvalue() <= 0.05 and 
                hinges[j].get_pvalue() <= 0.05 and 
                len(set(h1).intersection(h2)) > 0):
                hinge_union = list(set(h1).union(h2))
                
                flag = True
                for k in range(0, len(pre_overlap)):
                    if len(set(pre_overlap[k]).intersection(hinge_union)) > 0:
                        pre_overlap[k] = list(set(pre_overlap[k]).union(hinge_union))
                        flag = False
                        break
                
                if flag:
                    pre_overlap.append(hinge_union)
    
    # Merge overlapping hinges
    overlap = []
    for i in range(0, len(pre_overlap)):
        flag = True
        for j in range(0, len(overlap)):
            if len(set(overlap[j]).intersection(pre_overlap[i])) > 0:
                overlap[j] = list(set(overlap[j]).union(pre_overlap[i]))
                flag = False
        
        if flag:
            overlap.append(pre_overlap[i])
    
    # Filter to remove hinges at the end (within 3 residues)
    max_minus_3 = max([i.get_parent().get_id() for i in backbone_atoms]) - 3
    
    # Sort and filter
    overlap = sorted(overlap, key=lambda x: x[0])
    overlap = [sorted(i) for i in overlap if max_minus_3 not in i]
    
    return overlap, hinges

def page_home():
    """Home page"""
    st.title("🧬 PACKMAN Web Interface")
    st.write(f"**Version:** {get_version()}")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.header("About PACKMAN")
        st.write("""
        PACKMAN (Protein pACKing and Motion ANalysis) is a multi-utility tool 
        for studying protein packing and its effects on protein dynamics.
        
        ### Key Features:
        - 🔍 **Hinge Prediction**: Identify flexible protein hinges (domains)
        - 📊 **hdANM Analysis**: Elastic Network Model for dynamics
        - 📈 **Packing Entropy**: Voronoi-based entropy calculations
        """)
    
    with col2:
        st.header("Quick Links")
        st.markdown("""
        - 📚 [Documentation](https://py-packman.readthedocs.io/)
        - 🐙 [GitHub Repository](https://github.com/Pranavkhade/PACKMAN)
        """)
    
    st.divider()
    st.info("""
    **Getting Started:** Select a tool from the sidebar to begin your analysis.
    You can upload PDB files or use PDB IDs to fetch structures automatically.
    """)

def page_hinge_prediction():
    """Hinge Prediction page with ScanAlpha"""
    st.title("🔍 Hinge Prediction (ScanAlpha)")
    
    st.write("""
    This tool identifies protein hinges (flexible regions separating domains) 
    using an enhanced alpha scanning algorithm.
    """)
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("Input Settings")
        
        # File upload or PDB ID
        input_type = st.radio("Select input method:", ["Upload PDB", "PDB ID"])
        
        if input_type == "Upload PDB":
            uploaded_file = st.file_uploader("Choose a PDB/CIF file", type=["pdb", "cif"])
            if uploaded_file:
                pdb_path = uploaded_file.name
                file_extension = os.path.splitext(uploaded_file.name)[1].lower()
                with tempfile.NamedTemporaryFile(delete=False, suffix=f"{file_extension}") as f:
                    f.write(uploaded_file.read())
                    pdb_path = f.name
        else:
            pdb_id = st.text_input("Enter PDB ID (e.g., 1EXR):", value="1EXR")
            pdb_path = pdb_id
        
    
    with col2:
        st.subheader("Alpha Scan Parameters")
        alpha_start = st.number_input("Alpha Start:", min_value=0.0, max_value=1000.0, value=2.8, step=0.1)
        alpha_end = st.number_input("Alpha End:", min_value=0.0, max_value=10000.0, value=5.0, step=0.1)
        alpha_step = st.number_input("Alpha Step:", min_value=0.01, max_value=1.0, value=0.5, step=0.1)
    
    chain_id = st.text_input("Chain ID (leave empty for all chains):", value="")
    
    if st.button("Run Hinge Prediction", type="primary", use_container_width=True):
        try:
            with st.spinner("Loading structure..."):
                try:
                    mol = molecule.load_structure(pdb_path)
                except:
                    st.info("Downloading structure from PDB...")
                    try:
                        if input_type == "PDB ID":
                            molecule.download_structure(pdb_path, ftype='cif')
                            mol = molecule.load_structure(f"{pdb_path}.cif")
                        else:
                            mol = molecule.load_structure(pdb_path + '.cif')
                    except:
                        try:
                            molecule.download_structure(pdb_path, ftype='pdb')
                            mol = molecule.load_structure(f"{pdb_path}.pdb")
                        except Exception as e:
                            st.error(f"Could not load structure: {str(e)}")
                            return
            
            st.success("Structure loaded successfully!")
            
            # Get all hinges data
            all_hinges_data = []
            hinge_residues = []
            
            if chain_id:
                try:
                    backbone = [j for i in mol[0][chain_id].get_backbone() for j in i if j is not None]
                    overlap, hinges = scan_alpha_hinge_prediction(backbone, alpha_start, alpha_end, alpha_step)
                    
                    for idx, hinge_res in enumerate(overlap):
                        hinge_hinges = [h for h in hinges if 
                                       hinge_res[0] in [r.get_id() for r in h.get_elements()]]
                        p_value = hinge_hinges[0].get_pvalue() if hinge_hinges else "N/A"
                        all_hinges_data.append({
                            'ID': idx + 1,
                            'Chain': chain_id,
                            'Start': hinge_res[0],
                            'End': hinge_res[-1],
                            'P-value': p_value,
                            'color': hinge_colors[idx % len(hinge_colors)]
                        })
                        hinge_residues.append((hinge_res[0], hinge_res[-1]))
                except Exception as e:
                    st.error(f"Error processing chain {chain_id}: {str(e)}")
                    return
            else:
                st.info("Processing all chains...")
                for ch in mol[0].get_chains():
                    try:
                        backbone = [j for i in mol[0][ch.get_id()].get_backbone() for j in i if j is not None]
                        overlap, hinges = scan_alpha_hinge_prediction(backbone, alpha_start, alpha_end, alpha_step)
                        
                        for idx, hinge_res in enumerate(overlap):
                            hinge_hinges = [h for h in hinges if 
                                           hinge_res[0] in [r.get_id() for r in h.get_elements()]]
                            p_value = hinge_hinges[0].get_pvalue() if hinge_hinges else "N/A"
                            all_hinges_data.append({
                                'ID': len(all_hinges_data) + 1,
                                'Chain': ch.get_id(),
                                'Start': hinge_res[0],
                                'End': hinge_res[-1],
                                'P-value': p_value,
                                'color': hinge_colors[len(all_hinges_data) % len(hinge_colors)]
                            })
                            hinge_residues.append((hinge_res[0], hinge_res[-1]))
                    except:
                        continue

            # Display results and save to session state
            if all_hinges_data:
                st.success(f"Found {len(all_hinges_data)} hinge(s)!")
                
                # Save to session state
                st.session_state['hinge_results'] = all_hinges_data
                st.session_state['hinge_residues'] = hinge_residues
                st.session_state['hinge_pdb_path'] = pdb_path
            else:
                st.warning("No hinges detected with the given parameters. Try adjusting the alpha values.")
        
        except Exception as e:
            st.error(f"Error during analysis: {str(e)}")
            logger.exception("Hinge prediction error")

    # Display results from session state (persists across reruns)
    if 'hinge_results' in st.session_state and st.session_state['hinge_results']:
        all_hinges_data = st.session_state['hinge_results']
        hinge_residues = st.session_state['hinge_residues']
        pdb_path = st.session_state['hinge_pdb_path']
        
        # Display table
        st.subheader("Detected Hinges")
        st.dataframe(all_hinges_data, width="stretch")

        # Visualization
        st.subheader("3D Structure Visualization")

        try:
            render_protein_with_hinges(pdb_path, hinge_residues, chain_id if chain_id else None)
        except Exception as e:
            st.warning(f"Could not render 3D visualization: {str(e)}")
        
        # Download results
        st.subheader("Download .HNG file. (Necessary for hdANM hinge-based analysis)")
        
        # Select hinge numbers to download in .hng file
        selected_hinge_ids = st.multiselect("Select hinges to download:", [h['ID'] for h in all_hinges_data], key="hinge_ids")

        # chain lengths for .hng file
        chain_lengths = {}
        
        if( os.path.splitext(pdb_path)[-1] == '' ):
            mol = molecule.load_structure(pdb_path + '.cif')
        else:
            mol = molecule.load_structure(pdb_path)
            
        for ch in mol[0].get_chains():
            try:
                chain_lengths[ch.get_id()] = max([r.get_id() for r in mol[0][ch.get_id()].get_residues()])
            except:
                continue

        hng_content = generate_hng_file(pdb_path, all_hinges_data, selected_hinge_ids, chain_lengths)
        st.download_button(
            label="Download .HNG file",
            data=hng_content,
            file_name=f"{Path(pdb_path).stem}.hng",
            mime="text/plain"
        )

def generate_hng_file(pdb_path, hinges_data, selected_ids, chain_lengths):
    """Generate .hng format file"""
    filename = Path(pdb_path).stem
    content = []
    
    prev_end = 0
    domain_count = 1
    for i, hinge in enumerate(hinges_data):
        if hinge['ID'] in selected_ids:
            chain = hinge['Chain']
            start = hinge['Start']
            end = hinge['End']
            
            # add domain information
            if(i==0):
                content.append(f"{filename}_{chain}\tD{domain_count}\t1:{start - 1}")
                domain_count += 1
            elif(hinge['Chain'] != hinges_data[i-1]['Chain']):
                content.append(f"{filename}_{chain}\tD{domain_count}\t1:{start - 1}")
                domain_count += 1
            else:
                # start from previous hinge end + 1
                content.append(f"{filename}_{chain}\tD{domain_count}\t{prev_end + 1}:{start - 1}")
                domain_count += 1

            content.append(f"{filename}_{chain}\tH{hinge['ID']}\t{start}:{end}")

            # add domain information if another hinge doesnt exist in the given chain
            next_chain = hinges_data[i + 1]['Chain'] if i + 1 < len(hinges_data) else None
            if next_chain != chain:
                content.append(f"{filename}_{chain}\tD{domain_count}\t{end + 1}:{chain_lengths[chain]}")
                domain_count += 1
                # display warning if less than 3 residues in the last domain
                if(chain_lengths[chain] - end < 3):
                    st.warning(f"Warning: Last domain in chain {chain} has less than 3 residues, which may affect hdANM analysis.")
            
            prev_end = end
        
    return "\n".join(content)

def page_hdanm():
    """hdANM Analysis page"""
    st.title("📊 hdANM Analysis")
    
    st.write("""
    hdANM (Hinge-Domain Anisotropic Network Model) is a comprehensive 
    elastic network model for protein dynamics analysis. This tool generates 
    all modes of motion for your protein structure.
    """)
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("Input Settings")
        
        # File upload or PDB ID
        input_type = st.radio("Select input method:", ["Upload PDB", "PDB ID"], key="hdanm_input")
        
        if input_type == "Upload PDB":
            uploaded_file = st.file_uploader("Choose a PDB/CIF file", type=["pdb", "cif"], key="hdanm_file")
            if uploaded_file:
                pdb_path = uploaded_file.name
                with tempfile.NamedTemporaryFile(delete=False, suffix=f".{input_type.lower()}") as f:
                    f.write(uploaded_file.read())
                    pdb_path = f.name
        else:
            pdb_id = st.text_input("Enter PDB ID (e.g., 1EXR):", value="1EXR", key="hdanm_pdb_id")
            pdb_path = pdb_id
        
        uploaded_hng_file = st.file_uploader("Upload .HNG file for hinge-based analysis", type=["hng"], key="hdanm_hng_file")
        if uploaded_hng_file:
            with tempfile.NamedTemporaryFile(delete=False, suffix=f".{input_type.lower()}") as f:
                f.write(uploaded_hng_file.read())
                hng_path = f.name


    
    with col2:
        st.subheader("ANM Parameters")
        cutoff = st.number_input("Cutoff Distance (Å):", min_value=1.0, max_value=30.0, value=15.0)
        power = st.number_input("Power of Distance:", min_value=-5, max_value=5, value=0)
        mass_type = st.selectbox("Mass Type:", ["Molecular Weight", "Unit", "Atomics Mass"])
    
    chain_id = st.text_input("Chain ID (leave empty for all):", value="", key="hdanm_chain")
    
    if st.button("Run hdANM Analysis", type="primary", width=True, use_container_width=True):
        try:
            with st.spinner("Loading structure..."):
                try:
                    mol = molecule.load_structure(pdb_path)
                except:
                    st.info("Downloading structure from PDB...")
                    try:
                        if input_type == "PDB ID":
                            molecule.download_structure(pdb_path, ftype='cif')
                            mol = molecule.load_structure(f"{pdb_path}.cif")
                        else:
                            mol = molecule.load_structure(pdb_path + '.cif')
                    except:
                        try:
                            molecule.download_structure(pdb_path)
                            mol = molecule.load_structure(f"{pdb_path}.pdb")
                        except Exception as e:
                            st.error(f"Could not load structure: {str(e)}")
                            return
            
            st.success("Structure loaded!")
            
            # Get C-alpha atoms
            calpha = []
            if chain_id:
                try:
                    calpha = [i for i in mol[0][chain_id].get_calpha() if i is not None]
                except Exception as e:
                    st.error(f"Error processing chain {chain_id}: {str(e)}")
                    return
            else:
                for chain in mol[0].get_chains():
                    calpha.extend([i for i in mol[0][chain.get_id()].get_calpha() if i is not None])
            
            if not calpha:
                st.error("No C-alpha atoms found!")
                return
            
            with st.spinner("Computing hdANM..."):
                # Run hdANM
                Model = hdANM(calpha, hng_path, dr=cutoff, power=power)
                
                mass_map = {
                    "Molecular Weight": "residue",
                    "Unit": "unit",
                    "Atomics Mass": "atom"
                }
                Model.calculate_hessian(mass_type=mass_map[mass_type])
                Model.calculate_decomposition()
                
            
            st.success("Analysis complete!")
            
            st.session_state['hdanm_model'] = Model
            st.session_state['hdanm_eigenvalues'] = Model.get_eigenvalues()
            st.session_state['hdanm_eigenvectors'] = Model.get_eigenvectors()
            Model.calculate_fluctuations()
            st.session_state['hdanm_fluctuations'] = Model.get_fluctuations()
        except Exception as e:
            st.error(f"Error during hdANM analysis: {str(e)}")
            logger.exception("hdANM error")

    if 'hdanm_model' in st.session_state:
        # Display eigenvalues
        st.subheader("Eigenvalue Spectrum")
        eigenvalues = st.session_state['hdanm_eigenvalues']
        eigenvalues_sorted = sorted(eigenvalues)[:50]  # Show first 50
        # get the real part
        eigenvalues_sorted = [ev.real for ev in eigenvalues_sorted]
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.line_chart(eigenvalues_sorted)
        with col2:
            st.metric("Number of Modes", len(eigenvalues))
            st.metric("First Non-rigid Mode", "6" if len(eigenvalues) >= 6 else len(eigenvalues))
        with col3:
            fluctuations = st.session_state['hdanm_fluctuations']
            fluctuations = st.session_state['hdanm_fluctuations']/ max(st.session_state['hdanm_fluctuations'])
            st.line_chart([ev.real for ev in fluctuations])

        # Download options
        st.subheader("Download Results")
        col1, col2 = st.columns(2)
        
        with col1:
            twod_array = list(st.session_state['hdanm_eigenvectors'])

            #2d array to csv string
            full_eigenvectors_csv = ""
            for row in twod_array:
                full_eigenvectors_csv += ",".join([str(val.real) for val in row]) + "\n"
            

            st.download_button(
                label="Download Eigenvectors (CSV)",
                data=full_eigenvectors_csv,
                file_name="eigenvectors.csv",
                mime="text/csv"
            )
        
        with col2:
            eigenvalues_csv = np.array2string(st.session_state['hdanm_eigenvalues'], separator=',')
            st.download_button(
                label="Download Eigenvalues (CSV)",
                data=eigenvalues_csv,
                file_name="eigenvalues.csv",
                mime="text/csv"
            )
        
        # Movie generation
        st.subheader("Generate Motion Movies")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            mode_num = st.number_input("Mode Number:", min_value=6, max_value=len(eigenvalues), value=6)
        with col2:
            n_frames = st.number_input("Number of Frames:", min_value=5, max_value=50, value=10)
            st.session_state['hdanm_n_frames'] = n_frames
        with col3:
            scale = st.number_input("Movie Scale:", min_value=0.1, max_value=100.0, value=10.0)
        with col4:
            projection = st.selectbox("Projection Method:", ["Curvilinear", "Linear"])
        
        ca_to_aa = st.checkbox("Project C-Alpha motions to all atoms", value=True)

        if(ca_to_aa==False):
            st.warning("Movie will not appear if the all atoms are not projected. Please download and visualize .cif file locally.")
        
        if st.button("Generate Movie"):
            st.session_state['hdanm_movie_generated'] = False
            try:
                with st.spinner("Generating movie..."):
                    projection_map = {"Curvilinear": "curvilinear", "Linear": "linear"}
                    st.session_state['hdanm_model'].calculate_movie(
                        mode_num,
                        n=n_frames,
                        scale=scale,
                        extrapolation=projection_map[projection],
                        ca_to_aa=ca_to_aa
                    )
                
                st.success(f"Movie generated for mode {mode_num}!")

                st.session_state['hdanm_movie_file'] = f"{mode_num}.cif"
                st.session_state['hdanm_movie_generated'] = True
            except Exception as e:
                st.error(f"Error generating movie: {str(e)}")

        if 'hdanm_movie_generated' in st.session_state and st.session_state['hdanm_movie_generated']:
                st.subheader("3D Structure Visualization")

                try:
                    # render the generated movie inside a box
                    render_generated_movie(st.session_state['hdanm_movie_file'])
                except Exception as e:
                    st.warning(f"Could not render 3D visualization: {str(e)}")
                
                if os.path.exists(st.session_state['hdanm_movie_file']):
                    with open(st.session_state['hdanm_movie_file'], 'r') as f:
                        movie_content = f.read()
                        st.session_state['hdanm_movie_content'] = movie_content
                    st.download_button(
                        label="Download Movie (CIF)",
                        data=st.session_state['hdanm_movie_content'],
                        file_name=st.session_state['hdanm_movie_file'],
                        mime="text/plain"
                    )

def page_packing_entropy():
    """Voronoi Packing Entropy page"""
    st.title("📈 Voronoi Packing Entropy")
    
    st.write("""
    Calculate packing entropy using Voronoi tessellation to analyze 
    how tightly atoms are packed in your protein structure.
    """)
    st.info("Please remove hydrogen atoms from your structure for accurate entropy calculations, as they can skew the results.")
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("Input Settings")
        
        # File upload or PDB ID
        input_type = st.radio("Select input method:", ["Upload PDB", "PDB ID"], key="entropy_input")
        
        if input_type == "Upload PDB":
            uploaded_file = st.file_uploader("Choose a PDB/CIF file", type=["pdb", "cif"], key="entropy_file")
            if uploaded_file:
                pdb_path = uploaded_file.name
                
                extension = os.path.splitext(uploaded_file.name)[1].lower()
                with tempfile.NamedTemporaryFile(delete=False, suffix=f"{extension.lower()}") as f:
                    f.write(uploaded_file.read())
                    pdb_path = f.name
                    
        else:
            pdb_id = st.text_input("Enter PDB ID (e.g., 1EXR):", value="1EXR", key="entropy_pdb_id")
            pdb_path = pdb_id
    
    with col2:
        st.subheader("Entropy Parameters")
        probe_size = st.number_input("Probe Size (Å):", min_value=0.1, max_value=5.0, value=1.4, step=0.1)
        onsphere_points = st.number_input("On-Sphere Points:", min_value=10, max_value=100, value=30, step=5)
    
    chain_ids_input = st.text_input("Chain IDs (comma-separated, leave empty for all):", value="")
    
    if st.button("Calculate Packing Entropy", type="primary", use_container_width=True):
        try:
            if( os.path.splitext(pdb_path)[-1] == '' ):
                pdb_path = pdb_path + '.cif'
            else:
                pdb_path = pdb_path
            
            with st.spinner("Loading structure..."):
                try:
                    
                    if(os.path.splitext(pdb_path)[-1] == '.cif'):
                        mol = molecule.load_cif(pdb_path)
                    elif(os.path.splitext(pdb_path)[-1] == '.pdb'):
                        mol = molecule.load_pdb(pdb_path)
                    print(os.path.splitext(pdb_path)[-1])
                except:
                    st.info("Downloading structure from PDB...")
                    try:
                        if input_type == "PDB ID":
                            molecule.download_structure(pdb_path, ftype='cif')
                            mol = molecule.load_structure(f"{pdb_path}.cif")
                        else:
                            mol = molecule.load_structure(pdb_path + '.cif')
                    except:
                        try:
                            molecule.download_structure(pdb_path)
                            mol = molecule.load_structure(f"{pdb_path}.pdb")
                        except Exception as e:
                            st.error(f"Could not load structure: {str(e)}")
                            return
            
            st.success("Structure loaded!")
            
            # Parse chain IDs
            chains = None
            if chain_ids_input.strip():
                chains = [c.strip() for c in chain_ids_input.split(",")]
            
            with st.spinner("Computing packing entropy..."):
                results = []
                
                for frame in mol:
                    try:
                        result = PackingEntropy(
                            frame.get_atoms(),
                            chains=chains,
                            probe_size=probe_size,
                            onspherepoints=onsphere_points
                        )
                        
                        for residue in frame.get_residues():
                            try:
                                entropy = residue.get_entropy('PackingEntropy')
                                results.append({
                                    'Frame': frame.get_id(),
                                    'Chain': residue.get_parent().get_id(),
                                    'Residue_ID': residue.get_id(),
                                    'Residue_Name': residue.get_name(),
                                    'Packing_Entropy': entropy
                                })
                            except:
                                pass
                        
                        # Add frame totals
                        results.append({
                            'Frame': frame.get_id(),
                            'Chain': 'TOTAL',
                            'Residue_ID': None,
                            'Residue_Name': 'Frame Total',
                            'Packing_Entropy': result.get_total_entropy()
                        })
                        
                        normalized = result.get_total_entropy() / (3 * len([a for a in frame.get_atoms()]) - 6)
                        results.append({
                            'Frame': frame.get_id(),
                            'Chain': 'NORMALIZED',
                            'Residue_ID': None,
                            'Residue_Name': 'Normalized (3N-6)',
                            'Packing_Entropy': normalized
                        })
                    except Exception as e:
                        st.warning(f"Error processing frame: {str(e)}")
                        continue
            
            st.success("Calculation complete!")
            
            # Display results
            if results:
                st.subheader("Packing Entropy Results")
                st.dataframe(results, use_container_width=True)
                
                # Download results
                csv_content = "Frame,Chain,Residue_ID,Residue_Name,Packing_Entropy\n"
                for r in results:
                    csv_content += f"{r['Frame']},{r['Chain']},{r['Residue_ID']},{r['Residue_Name']},{r['Packing_Entropy']}\n"
                
                
            else:
                st.warning("No results obtained. Check your parameters.")
        
        except Exception as e:
            st.error(f"Error during calculation: {str(e)}")
            logger.exception("Packing entropy error")

def main():
    """Main application"""
    # Sidebar navigation
    with st.sidebar:
        st.markdown("# PACKMAN Tools")
        st.divider()
        
        page = st.radio(
            "Select Analysis:",
            [
                "🏠 Home",
                "🔍 Hinge Prediction",
                "📊 hdANM Analysis",
                "📈 Packing Entropy"
            ]
        )
        
        st.divider()
        st.markdown("### About")
        st.info(f"PACKMAN v{get_version()}")
        st.markdown("""
        [GitHub](https://github.com/Pranavkhade/PACKMAN) • 
        [Docs](https://py-packman.readthedocs.io/)""")
    
    # Route to pages
    if page == "🏠 Home":
        page_home()
    elif page == "🔍 Hinge Prediction":
        page_hinge_prediction()
    elif page == "📊 hdANM Analysis":
        page_hdanm()
    elif page == "📈 Packing Entropy":
        page_packing_entropy()

if __name__ == "__main__":
    main()
