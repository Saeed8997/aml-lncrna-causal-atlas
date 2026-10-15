"""
src/multi_resolution.py
Multi-resolution causal networks: Global, LSC-like, and Monocytic blast states with differential rewiring.
"""

import pandas as pd
import numpy as np
from src.causal_discovery import run_bootstrap_pc_discovery
from src.utils import setup_logger

logger = setup_logger("multi_resolution")

def run_multi_resolution_causal_discovery(
    meta_adata,
    gene_types: dict,
    config: dict
) -> dict:
    """
    Executes causal discovery across three biological resolutions:
    1. Global AML network
    2. LSC-like / immature blast state
    3. Differentiated monocytic blast state
    """
    results = {}
    df_all = pd.DataFrame(meta_adata.X, columns=meta_adata.var_names)
    
    # 1. Global AML Causal Network
    logger.info("--- Resolution 1: Global AML Causal Network ---")
    global_edges = run_bootstrap_pc_discovery(df_all, gene_types, config)
    global_edges["network_resolution"] = "Global"
    results["global"] = global_edges
    
    # 2. LSC-like Immature Blast Network
    lsc_mask = (meta_adata.obs["macro_state"] == "LSC_like").values
    if np.sum(lsc_mask) >= 1:
        logger.info(f"--- Resolution 2: LSC-like State Network (N={np.sum(lsc_mask)} metacells) ---")
        df_lsc = df_all[lsc_mask].reset_index(drop=True)
        lsc_edges = run_bootstrap_pc_discovery(df_lsc, gene_types, config)
        lsc_edges["network_resolution"] = "LSC_like"
        results["lsc"] = lsc_edges
    else:
        logger.warning("Insufficient LSC-like metacells; using global as fallback.")
        results["lsc"] = global_edges.copy()
        results["lsc"]["network_resolution"] = "LSC_like"

    # 3. Monocytic Differentiated Blast Network
    mono_mask = (meta_adata.obs["macro_state"] == "Monocytic").values
    if np.sum(mono_mask) >= 1:
        logger.info(f"--- Resolution 3: Monocytic State Network (N={np.sum(mono_mask)} metacells) ---")
        df_mono = df_all[mono_mask].reset_index(drop=True)
        mono_edges = run_bootstrap_pc_discovery(df_mono, gene_types, config)
        mono_edges["network_resolution"] = "Monocytic"
        results["mono"] = mono_edges
    else:
        logger.warning("Insufficient Monocytic metacells; using global as fallback.")
        results["mono"] = global_edges.copy()
        results["mono"]["network_resolution"] = "Monocytic"

    return results

def compute_differential_causal_rewiring(
    lsc_edges: pd.DataFrame,
    mono_edges: pd.DataFrame,
    all_genes: list
) -> tuple:
    """
    Computes edge-level and gene-level rewiring between LSC-like and Monocytic states.
    Edge weight = bootstrap_freq * correlation.
    Delta W = |W_LSC - W_Mono|.
    """
    logger.info("Computing Differential Causal Rewiring between LSC and Monocytic states...")
    
    weights_lsc = {}
    for _, row in lsc_edges.iterrows():
        key = (row["source"], row["target"])
        weights_lsc[key] = row["bootstrap_freq"] * row["correlation"]
        
    weights_mono = {}
    for _, row in mono_edges.iterrows():
        key = (row["source"], row["target"])
        weights_mono[key] = row["bootstrap_freq"] * row["correlation"]
        
    all_edge_keys = set(weights_lsc.keys()) | set(weights_mono.keys())
    
    rewiring_rows = []
    gene_rewiring_sum = {g: 0.0 for g in all_genes}
    
    for u, v in all_edge_keys:
        w_l = weights_lsc.get((u, v), 0.0)
        w_m = weights_mono.get((u, v), 0.0)
        delta = abs(w_l - w_m)
        
        state_class = "Conserved"
        if (u, v) in weights_lsc and (u, v) not in weights_mono:
            state_class = "LSC_Specific"
        elif (u, v) in weights_mono and (u, v) not in weights_lsc:
            state_class = "Monocytic_Specific"
        elif np.sign(w_l) != np.sign(w_m):
            state_class = "Sign_Inversion"
            
        rewiring_rows.append({
            "source": u,
            "target": v,
            "weight_lsc": float(w_l),
            "weight_mono": float(w_m),
            "delta_weight": float(delta),
            "rewiring_class": state_class
        })
        
        gene_rewiring_sum[u] = gene_rewiring_sum.get(u, 0.0) + delta
        
    rewiring_df = pd.DataFrame(rewiring_rows)
    if not rewiring_df.empty:
        rewiring_df = rewiring_df.sort_values(by="delta_weight", ascending=False).reset_index(drop=True)
        
    gene_scores = pd.DataFrame([
        {"gene": g, "total_rewiring_magnitude": gene_rewiring_sum.get(g, 0.0)}
        for g in all_genes
    ]).sort_values(by="total_rewiring_magnitude", ascending=False).reset_index(drop=True)
    
    logger.info(f"Identified {len(rewiring_df)} differential edges across {len(gene_scores)} genes.")
    return rewiring_df, gene_scores
