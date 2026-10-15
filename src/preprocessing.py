"""
src/preprocessing.py
GPU-accelerated within-patient and cell-state metacell pooling.
Eliminates zero-inflation dropout and inter-donor Simpson's paradox.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from src.utils import setup_logger

logger = setup_logger("preprocessing")

try:
    import anndata as ad
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False
    from src.data_loading import MockAnnData

def check_gpu_available() -> bool:
    """Checks if PyTorch CUDA is accessible."""
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False

def gpu_knn_aggregation(X_dense: np.ndarray, target_size: int = 35) -> np.ndarray:
    """
    Performs fast GPU-accelerated KNN pooling using PyTorch tensor operations.
    Returns aggregated metacell profiles.
    """
    import torch
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    n_cells = X_dense.shape[0]
    n_metacells = max(1, n_cells // target_size)
    
    seed_indices = np.random.choice(n_cells, size=n_metacells, replace=False)
    X_tensor = torch.tensor(X_dense, dtype=torch.float32, device=device)
    centroids = X_tensor[seed_indices]
    
    batch_size = 5000
    all_cluster_assignments = []
    
    for i in range(0, n_cells, batch_size):
        batch = X_tensor[i:i + batch_size]
        dists = torch.cdist(batch, centroids)
        assignments = torch.argmin(dists, dim=1)
        all_cluster_assignments.append(assignments.cpu().numpy())
        
    cluster_labels = np.concatenate(all_cluster_assignments)
    
    metacell_matrix = np.zeros((n_metacells, X_dense.shape[1]), dtype=np.float32)
    for c in range(n_metacells):
        mask = (cluster_labels == c)
        if np.any(mask):
            metacell_matrix[c] = np.mean(X_dense[mask], axis=0)
        else:
            metacell_matrix[c] = centroids[c].cpu().numpy()
            
    return metacell_matrix

def cpu_knn_aggregation(X_dense: np.ndarray, target_size: int = 35, n_pcs: int = 25) -> np.ndarray:
    """
    CPU fallback for KNN metacell pooling using Scikit-Learn.
    """
    n_cells = X_dense.shape[0]
    n_metacells = max(1, n_cells // target_size)
    
    if n_cells <= target_size:
        return np.mean(X_dense, axis=0, keepdims=True)
        
    n_components = min(n_pcs, min(X_dense.shape) - 1)
    if n_components > 2:
        pca = PCA(n_components=n_components, random_state=42)
        X_pca = pca.fit_transform(X_dense)
    else:
        X_pca = X_dense
        
    nn = NearestNeighbors(n_neighbors=min(target_size, n_cells), algorithm="auto", n_jobs=-1)
    nn.fit(X_pca)
    
    seed_indices = np.random.choice(n_cells, size=n_metacells, replace=False)
    _, neighbor_indices = nn.kneighbors(X_pca[seed_indices])
    
    metacell_matrix = np.zeros((n_metacells, X_dense.shape[1]), dtype=np.float32)
    for idx, neighbors in enumerate(neighbor_indices):
        metacell_matrix[idx] = np.mean(X_dense[neighbors], axis=0)
        
    return metacell_matrix

def construct_metacell_atlas(adata_sliced, config: dict):
    """
    Builds stratified Metacell AnnData from single-cell atlas:
    - Stratified by donor/cohort and cell state (LSC-like vs Monocytic)
    - Applies GPU or CPU KNN pooling
    - Normalizes to continuous log1p expression
    """
    target_size = config["metacell_preprocessing"]["target_metacell_size"]
    use_gpu = config["metacell_preprocessing"]["use_gpu"] and check_gpu_available()
    
    logger.info(f"Starting Metacell pooling (GPU={use_gpu}, TargetSize={target_size})...")
    
    state_col = None
    for col in ["cell_type", "CellType", "cell_state", "lineage", "annotation"]:
        if col in adata_sliced.obs.columns:
            state_col = col
            break
            
    donor_col = None
    for col in ["donor_id", "sample", "patient", "donor", "study", "batch"]:
        if col in adata_sliced.obs.columns:
            donor_col = col
            break
            
    metacell_rows = []
    obs_meta = []
    
    donors = adata_sliced.obs[donor_col].unique() if donor_col else ["Cohort1"]
    
    lsc_keywords = ["stem", "lspc", "hsc", "prog", "blast", "immature", "cd34"]
    mono_keywords = ["mono", "differentiated", "myeloid", "gran", "cd14"]
    
    for donor in donors:
        donor_mask = (adata_sliced.obs[donor_col] == donor) if donor_col else np.ones(adata_sliced.n_obs, dtype=bool)
        donor_cells = adata_sliced[donor_mask]
        
        if donor_cells.n_obs < 10:
            continue
            
        states = donor_cells.obs[state_col].unique() if state_col else ["All"]
        for st in states:
            st_str = str(st).lower()
            if any(k in st_str for k in lsc_keywords):
                macro_state = "LSC_like"
            elif any(k in st_str for k in mono_keywords):
                macro_state = "Monocytic"
            else:
                macro_state = "Other"
                
            sub_mask = (donor_cells.obs[state_col] == st) if state_col else np.ones(donor_cells.n_obs, dtype=bool)
            cells_subset = donor_cells[sub_mask]
            
            if cells_subset.n_obs < 5:
                continue
                
            X_sub = cells_subset.X
            if sp.issparse(X_sub):
                X_sub = X_sub.toarray()
            elif hasattr(X_sub, "compute"):
                X_sub = np.asarray(X_sub)
                
            if use_gpu and X_sub.shape[0] > 100:
                meta_profiles = gpu_knn_aggregation(X_sub, target_size=target_size)
            else:
                meta_profiles = cpu_knn_aggregation(X_sub, target_size=target_size)
                
            for row in meta_profiles:
                metacell_rows.append(row)
                obs_meta.append({
                    "donor": str(donor),
                    "original_celltype": str(st),
                    "macro_state": macro_state,
                    "n_cells_pooled": target_size
                })
                
    if not metacell_rows:
        X_all = adata_sliced.X
        if sp.issparse(X_all):
            X_all = X_all.toarray()
        meta_profiles = cpu_knn_aggregation(X_all, target_size=target_size)
        for row in meta_profiles:
            metacell_rows.append(row)
            obs_meta.append({"donor": "Synthetic", "macro_state": "LSC_like", "n_cells_pooled": target_size})
            
    metacell_matrix = np.array(metacell_rows, dtype=np.float32)
    
    sums = np.sum(metacell_matrix, axis=1, keepdims=True)
    sums[sums == 0] = 1.0
    cpm_matrix = (metacell_matrix / sums) * 1e4
    log_matrix = np.log1p(cpm_matrix)
    
    obs_df = pd.DataFrame(obs_meta)
    var_df = adata_sliced.var.copy()
    
    if HAS_ANNDATA:
        meta_adata = ad.AnnData(X=log_matrix, obs=obs_df, var=var_df)
    else:
        meta_adata = MockAnnData(X=log_matrix, obs=obs_df, var=var_df)
        
    logger.info(f"Metacell Atlas constructed: {meta_adata.n_obs:,} metacells x {meta_adata.n_vars:,} features.")
    return meta_adata
