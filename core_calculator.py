import numpy as np
import scipy.spatial.distance as dist

class CytoSphereEngine:
    def __init__(self, ligand_expr, receptor_expr, spatial_coords, morphology_scores):
        """
        Initializes the engine with 1D arrays of expression, 2D coords, and 1D morpho scores.
        """
        self.L = np.array(ligand_expr)
        self.R = np.array(receptor_expr)
        self.coords = np.array(spatial_coords)
        self.M = np.array(morphology_scores)
        self.n_cells = len(self.L)

    def calculate_spatial_decay(self, beta=50.0, max_radius=150.0):
        """
        Calculates delta(d_i,j) = e^(-d / beta)
        Cuts off signaling entirely if distance > max_radius.
        """
        # Calculate pairwise Euclidean distances
        distances = dist.cdist(self.coords, self.coords, 'euclidean')
        
        # Apply exponential decay
        decay_matrix = np.exp(-distances / beta)
        
        # Zero out distances beyond biologically feasible paracrine signaling
        decay_matrix[distances > max_radius] = 0
        
        # Remove autocrine signaling (diagonal = 0)
        np.fill_diagonal(decay_matrix, 0)
        return decay_matrix

    def compute_cip(self, beta=50.0, max_radius=150.0):
        """
        Executes the novel CIP formula.
        """
        # 1. Transcriptomic baseline: L_i * R_j (Outer product)
        transcriptomic_base = np.outer(self.L, self.R)
        
        # 2. Apply Spatial Decay: delta(d)
        spatial_decay_matrix = self.calculate_spatial_decay(beta, max_radius)
        spatio_transcriptomic = transcriptomic_base * spatial_decay_matrix
        
        # 3. Apply Morphological Mask: mu_j
        # Tile the morphology scores to multiply against the receiver columns
        morpho_matrix = np.tile(self.M, (self.n_cells, 1))
        
        # Final CytoSphere Interaction Probability Matrix
        cip_matrix = spatio_transcriptomic * morpho_matrix
        
        return cip_matrix
      
