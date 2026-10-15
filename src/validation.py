"""
src/validation.py
Statistical validation against GSE318789 functional CRISPRi screens (proliferation, differentiation, Venetoclax).
"""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom
from src.utils import setup_logger

logger = setup_logger("validation")

def run_crispri_enrichment_validation(
    mrs_df: pd.DataFrame,
    crispri_hits_dict: dict,
    top_k_threshold: int = 25
) -> pd.DataFrame:
    """
    Performs Fisher's Exact Test and hypergeometric enrichment for top causal master regulators.
    Evaluates across Proliferation, Differentiation, and Venetoclax sensitivity screens.
    """
    logger.info("Executing Statistical Validation against GSE318789 CRISPRi screens...")
    
    # Filter to lncRNAs in network
    lncrna_df = mrs_df[mrs_df["gene_type"] == "lncrna"].copy().reset_index(drop=True)
    N = len(lncrna_df)
    if N == 0:
        logger.warning("No lncRNAs identified in MRS table; evaluating all genes.")
        lncrna_df = mrs_df.copy()
        N = len(lncrna_df)
        
    top_k = min(top_k_threshold, N)
    top_genes = set(lncrna_df.iloc[:top_k]["gene"])
    non_top_genes = set(lncrna_df.iloc[top_k:]["gene"])
    
    validation_records = []
    
    for screen_name, hit_list in crispri_hits_dict.items():
        hits_in_network = set(hit_list) & set(lncrna_df["gene"])
        K = len(hits_in_network)
        
        # 2x2 contingency table:
        #                 Top_MRS    Non_Top_MRS
        # CRISPRi_Hit        a            b
        # Not_Hit            c            d
        
        a = len(top_genes & hits_in_network)
        b = len(non_top_genes & hits_in_network)
        c = len(top_genes - hits_in_network)
        d = len(non_top_genes - hits_in_network)
        
        # Contingency table
        table = [[a, b], [c, d]]
        odds_ratio, p_val = fisher_exact(table, alternative="greater")
        
        # Hypergeometric test: P(X >= a)
        # N: population, K: success states in population, n: draws (top_k)
        hyper_p = hypergeom.sf(a - 1, N, K, top_k) if K > 0 else 1.0
        
        expected = (K / N) * top_k if N > 0 else 0
        fold_enrichment = (a / expected) if expected > 0 else 0.0
        
        overlapping_hits = list(top_genes & hits_in_network)
        
        validation_records.append({
            "screen_phenotype": screen_name,
            "total_network_lncrnas": N,
            "crispri_hits_in_network": K,
            "top_k_evaluated": top_k,
            "top_k_crispri_hits": a,
            "expected_hits": round(expected, 2),
            "fold_enrichment": round(fold_enrichment, 2),
            "odds_ratio": round(odds_ratio, 2) if not np.isinf(odds_ratio) else 999.0,
            "fisher_pvalue": float(p_val),
            "hypergeom_pvalue": float(hyper_p),
            "validated_hit_symbols": ", ".join(overlapping_hits) if overlapping_hits else "None"
        })
        
        logger.info(
            f"Screen '{screen_name}': {a}/{top_k} hits (Fold={fold_enrichment:.2f}x, "
            f"Fisher p={p_val:.4e}, OR={odds_ratio:.2f}) -> Overlap: {overlapping_hits}"
        )
        
    val_df = pd.DataFrame(validation_records)
    return val_df
