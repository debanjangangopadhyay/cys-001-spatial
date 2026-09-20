import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from core_calculator import CytoSphereEngine

st.set_page_config(page_title="CytoSphere AI", layout="wide")

st.title("CytoSphere AI: Spatio-Morphological Communication Inference")
st.markdown("A novel framework integrating spatial transcriptomics with cytological morphology.")

# Sidebar Controls
st.sidebar.header("Parameters")
beta_val = st.sidebar.slider("Signal Decay Factor (β)", 10.0, 100.0, 50.0)
threshold = st.sidebar.slider("Interaction Threshold", 0.0, 1.0, 0.5)

st.write("### Simulation Mode")
st.info("For demonstration, we are generating synthetic cell data representing an H&E slide coordinate system.")

# Generate Synthetic Data for the App
if st.button("Generate & Analyze Tissue Map"):
    n_cells = 50
    np.random.seed(42)
    
    # Synthetic Data: Coordinates, Expression, Atypia
    coords = np.random.rand(n_cells, 2) * 500  # 500x500 pixel area
    ligand = np.random.rand(n_cells)
    receptor = np.random.rand(n_cells)
    morphology = np.random.rand(n_cells) # E.g., Nuclear Atypia Score
    
    # Initialize Engine
    engine = CytoSphereEngine(ligand, receptor, coords, morphology)
    cip_matrix = engine.compute_cip(beta=beta_val)
    
    # Plotting the Novel Network
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Draw cells
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=morphology, cmap='Reds', s=100, edgecolors='black')
    plt.colorbar(scatter, label="Morphological Atypia Score (μ)")
    
    # Draw interactions above threshold
    sources, targets = np.where(cip_matrix > threshold)
    for src, tgt in zip(sources, targets):
        x = [coords[src, 0], coords[tgt, 0]]
        y = [coords[src, 1], coords[tgt, 1]]
        weight = cip_matrix[src, tgt]
        ax.plot(x, y, color='cyan', linewidth=weight*5, alpha=0.6)
        
    ax.set_title("In-Situ CytoSphere Network Overlay")
    ax.set_facecolor('#222222') # Dark background to simulate fluorescence/slide contrast
    st.pyplot(fig)
    
    st.success("Calculated CytoSphere Interaction Probabilities successfully.")

