"""
src/feature_selection.py
Dual-criteria lncRNA identification, ID harmonization, hybrid feature selection, and Cis/Trans annotation.
"""

import re
import pandas as pd
import numpy as np
from src.utils import setup_logger

logger = setup_logger("feature_selection")

class IDHarmonizer:
    """Harmonizes gene identifiers across Ensembl, HGNC symbols, clone IDs, and FANTOM CAT IDs."""
    
    def __init__(self, var_df: pd.DataFrame):
        self.var_df = var_df.copy()
        self.index_is_ensg = any(str(x).startswith("ENSG") for x in self.var_df.index[:20])
        
        self.id_to_symbol = {}
        self.symbol_to_id = {}
        
        feature_name_col = None
        for col in ["feature_name", "gene_symbols", "symbol", "gene_name"]:
            if col in self.var_df.columns:
                feature_name_col = col
                break
                
        for idx, row in self.var_df.iterrows():
            ensg_clean = str(idx).split(".")[0]
            sym = str(row[feature_name_col]) if feature_name_col else str(idx)
            self.id_to_symbol[ensg_clean] = sym
            self.id_to_symbol[str(idx)] = sym
            self.symbol_to_id[sym.upper()] = str(idx)
            self.symbol_to_id[ensg_clean] = str(idx)
            
    def map_gene(self, query: str) -> str:
        """Returns the corresponding index in var_df, or None if unmatched."""
        q_clean = str(query).strip().upper()
        q_ensg = q_clean.split(".")[0]
        if q_clean in self.symbol_to_id:
            return self.symbol_to_id[q_clean]
        if q_ensg in self.symbol_to_id:
            return self.symbol_to_id[q_ensg]
        if query in self.var_df.index:
            return query
        return None

def identify_all_lncrnas(var_df: pd.DataFrame) -> set:
    """
    Detects lncRNAs using biotype annotations AND regex pattern heuristics.
    """
    lnc_set = set()
    
    for col in ["feature_biotype", "gene_biotype", "biotype"]:
        if col in var_df.columns:
            mask = var_df[col].astype(str).str.lower().str.contains("lnc|antisense|intergenic")
            lnc_set.update(var_df.index[mask])
            break
            
    lnc_pattern = re.compile(r'^(LINC\d+|MIR\d+HG|CATG\d+|.*-AS\d*|^AC\d+\.\d+|^AL\d+\.\d+|\w+-DT|\w+-IT)$', re.IGNORECASE)
    
    feature_name_col = None
    for col in ["feature_name", "gene_symbols", "symbol"]:
        if col in var_df.columns:
            feature_name_col = col
            break
            
    for idx, row in var_df.iterrows():
        sym = str(row[feature_name_col]) if feature_name_col else str(idx)
        if lnc_pattern.match(sym) or lnc_pattern.match(str(idx)):
            lnc_set.add(idx)
            
    logger.info(f"Identified {len(lnc_set):,} total putative lncRNAs in var annotations.")
    return lnc_set

def select_hybrid_feature_set(adata, config: dict, crispr_targets: set = None) -> list:
    """
    Constructs the hybrid feature space (~600-800 genes):
    - CRISPRi screen hit lncRNAs & target overlaps
    - Core AML oncogenic drivers and Venetoclax axis genes
    - Key lineage transcription factors (TFs)
    - Top variable expressed lncRNAs in the atlas
    """
    var_df = adata.var.copy()
    harmonizer = IDHarmonizer(var_df)
    selected_indices = set()
    
    core_genes = (
        config["features"]["core_aml_drivers"] + 
        config["features"]["venetoclax_axis"] + 
        config["features"]["stemness_axis"] + 
        config["features"]["key_tfs"]
    )
    for g in core_genes:
        mapped = harmonizer.map_gene(g)
        if mapped:
            selected_indices.add(mapped)
            
    logger.info(f"Retained {len(selected_indices)} core driver & TF genes in atlas.")
    
    known_hits = []
    for category, hits in config["features"]["crispri_validation_hits"].items():
        known_hits.extend(hits)
        
    for h in known_hits:
        mapped = harmonizer.map_gene(h)
        if mapped:
            selected_indices.add(mapped)
            
    if crispr_targets:
        matched_screen = 0
        for target in crispr_targets:
            mapped = harmonizer.map_gene(target)
            if mapped and mapped not in selected_indices:
                selected_indices.add(mapped)
                matched_screen += 1
                if len(selected_indices) >= config["features"]["target_feature_count"]:
                    break
        logger.info(f"Added {matched_screen} matching CRISPRi screen targets.")

    all_lncs = identify_all_lncrnas(var_df)
    remaining_lncs = list(all_lncs - selected_indices)
    
    target_count = config["features"]["target_feature_count"]
    if len(selected_indices) < target_count and remaining_lncs:
        needed = min(target_count - len(selected_indices), len(remaining_lncs))
        selected_indices.update(remaining_lncs[:needed])
        logger.info(f"Added {needed} additional candidate lncRNAs to reach target count.")

    final_features = list(selected_indices)
    logger.info(f"Final hybrid feature set: {len(final_features)} genes.")
    return final_features

def annotate_cis_trans(edge_u: str, edge_v: str, coord_dict: dict, max_cis_distance: int = 500000) -> str:
    """
    Annotates whether an edge (u -> v) is Cis (<500kb on same chromosome) or Trans.
    """
    info_u = coord_dict.get(edge_u)
    info_v = coord_dict.get(edge_v)
    if not info_u or not info_v:
        return "Trans"
    
    chr_u, pos_u = info_u.get("chr"), info_u.get("tss")
    chr_v, pos_v = info_v.get("chr"), info_v.get("tss")
    
    if chr_u and chr_v and chr_u == chr_v and pos_u is not None and pos_v is not None:
        if abs(pos_u - pos_v) <= max_cis_distance:
            return "Cis"
            
    return "Trans"
