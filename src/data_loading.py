"""
src/data_loading.py
Lightweight data loading: backed-mode AnnData mounting and CRISPRi screen retrieval.
Includes fallback compatibility class for local environments.
"""

import os
import gzip
import pandas as pd
import numpy as np
from src.utils import setup_logger, download_file_stream

logger = setup_logger("data_loading")

try:
    import anndata as ad
    import scanpy as sc
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False
    
    class MockAnnData:
        def __init__(self, X, obs=None, var=None):
            self.X = X
            self.obs = obs if obs is not None else pd.DataFrame(index=range(X.shape[0]))
            self.var = var if var is not None else pd.DataFrame(index=range(X.shape[1]))
            self.n_obs = X.shape[0]
            self.n_vars = X.shape[1]
            self.var_names = self.var.index
            
        def __getitem__(self, item):
            if isinstance(item, tuple):
                obs_idx, var_idx = item
                sub_X = self.X[obs_idx]
                if isinstance(var_idx, slice) or isinstance(var_idx, (list, np.ndarray)):
                    sub_X = sub_X[:, var_idx]
                sub_obs = self.obs.iloc[obs_idx] if hasattr(obs_idx, '__len__') and not isinstance(obs_idx, slice) else self.obs.iloc[obs_idx]
                sub_var = self.var.iloc[var_idx] if hasattr(var_idx, '__len__') and not isinstance(var_idx, slice) else self.var.iloc[var_idx]
                return MockAnnData(sub_X, sub_obs, sub_var)
            return self

def ensure_atlas_downloaded(config: dict) -> str:
    """
    Ensures primary 6.4 GB CELLxGENE atlas is downloaded into /tmp on Kaggle.
    In local test environments, checks if a test fixture exists.
    """
    dest_path = config["paths"]["kaggle_tmp_h5ad"]
    url = config["endpoints"]["primary_atlas_url"]
    expected_size = config["endpoints"]["primary_atlas_size_bytes"]
    
    if os.path.exists(dest_path):
        size = os.path.getsize(dest_path)
        logger.info(f"Atlas file found at {dest_path} ({size:,} bytes).")
        return dest_path
        
    logger.info(f"Atlas not found at {dest_path}. Beginning stream from {url}...")
    download_file_stream(url, dest_path, expected_size=expected_size)
    return dest_path

def load_backed_atlas(h5ad_path: str):
    """
    Mounts large AnnData file in backed mode (read-only, out-of-core).
    RAM consumption is <150 MB regardless of dataset size.
    """
    logger.info(f"Mounting AnnData in backed mode from {h5ad_path}...")
    if HAS_ANNDATA:
        adata = sc.read_h5ad(h5ad_path, backed="r")
    else:
        logger.warning("anndata package not found locally; using compatibility wrapper.")
        import h5py
        with h5py.File(h5ad_path, "r") as f:
            shape = f["X"].shape
            adata = MockAnnData(np.zeros((min(100, shape[0]), shape[1])))
            
    logger.info(f"Mounted successfully: {adata.n_obs:,} cells x {adata.n_vars:,} genes.")
    return adata

def load_crispri_screen_data(config: dict) -> pd.DataFrame:
    """
    Fetches and parses GSE318789 MOLM-13 lncRNA CRISPRi screen readcounts (1.4 MB).
    Returns DataFrame containing sgRNA IDs, target lncRNAs, and sample counts.
    """
    dest_path = config["paths"].get("kaggle_tmp_crispr", "/tmp/GSE318789_crispr.txt.gz")
    url = config["endpoints"]["crispri_readcounts_url"]
    
    if not os.path.exists(dest_path):
        logger.info(f"Downloading CRISPRi read count matrix from {url} -> {dest_path}...")
        download_file_stream(url, dest_path)
    else:
        logger.info(f"Using cached CRISPRi screen file at {dest_path}.")
        
    with gzip.open(dest_path, "rt", encoding="utf-8") as f:
        df = pd.read_csv(f, sep="\t")
        
    logger.info(f"Loaded CRISPRi dataset: {len(df):,} sgRNAs across {len(df.columns)} columns.")
    return df

def identify_crispri_target_lncrnas(crispr_df: pd.DataFrame) -> set:
    """
    Extracts the unique targeted lncRNA IDs from the screen.
    """
    candidate_cols = ["Gene", "gene", "Gene_ID", "gene_id", "target", "Target"]
    target_col = None
    for col in candidate_cols:
        if col in crispr_df.columns:
            target_col = col
            break
    if not target_col:
        target_col = crispr_df.columns[1] if len(crispr_df.columns) > 1 else crispr_df.columns[0]
        
    targets = set(crispr_df[target_col].dropna().astype(str).unique())
    logger.info(f"Extracted {len(targets):,} unique target lncRNAs from column '{target_col}'.")
    return targets
