"""
src/graph_analysis.py
Graph centrality, Causal Cascade Reach, Master Regulator Score (MRS), and GraphML export.
"""

import numpy as np
import pandas as pd
import networkx as nx
from src.utils import setup_logger

logger = setup_logger("graph_analysis")

def build_networkx_graph(edges_df: pd.DataFrame, directed_only: bool = True) -> nx.DiGraph:
    """Builds a NetworkX DiGraph from stable edge table."""
    G = nx.DiGraph()
    for _, row in edges_df.iterrows():
        if directed_only and not row.get("directed", True):
            continue
        u, v = row["source"], row["target"]
        G.add_edge(
            u, v,
            weight=float(row.get("bootstrap_freq", 1.0)),
            regulation_sign=int(row.get("regulation_sign", 1)),
            dir_prob=float(row.get("dir_prob", 1.0)),
            cis_trans=str(row.get("cis_trans", "Trans"))
        )
    return G

def compute_master_regulator_scores(
    stable_edges: pd.DataFrame,
    gene_scores: pd.DataFrame,
    gene_types: dict
) -> pd.DataFrame:
    """
    Computes composite Master Regulator Score (MRS):
    MRS = 0.35*OutDegree_norm + 0.30*CascadeReach_norm + 0.20*Stability_mean + 0.15*Rewiring_norm
    """
    logger.info("Computing Master Regulator Scores (MRS) across gene regulatory networks...")
    
    G = build_networkx_graph(stable_edges, directed_only=True)
    all_nodes = set(stable_edges["source"]) | set(stable_edges["target"])
    for n in all_nodes:
        if n not in G:
            G.add_node(n)
            
    rewiring_map = dict(zip(gene_scores["gene"], gene_scores["total_rewiring_magnitude"])) if gene_scores is not None else {}
    
    mrs_records = []
    for node in all_nodes:
        out_deg = G.out_degree(node)
        in_deg = G.in_degree(node)
        
        # Causal Cascade Reach (downstream transitive descendants)
        try:
            descendants = nx.descendants(G, node)
            cascade_reach = len(descendants)
        except Exception:
            cascade_reach = out_deg
            
        # Outgoing edge stability mean
        out_edges = stable_edges[(stable_edges["source"] == node) & (stable_edges["directed"])]
        stability_mean = float(out_edges["bootstrap_freq"].mean()) if not out_edges.empty else 0.0
        
        rewiring_val = float(rewiring_map.get(node, 0.0))
        g_type = gene_types.get(node, "effector")
        
        mrs_records.append({
            "gene": node,
            "gene_type": g_type,
            "out_degree": out_deg,
            "in_degree": in_deg,
            "cascade_reach": cascade_reach,
            "stability_mean": stability_mean,
            "rewiring_magnitude": rewiring_val
        })
        
    df_mrs = pd.DataFrame(mrs_records)
    
    # Normalize features to [0, 1]
    def _minmax(s):
        denom = s.max() - s.min()
        return (s - s.min()) / denom if denom > 0 else np.zeros_like(s)
        
    df_mrs["norm_out"] = _minmax(df_mrs["out_degree"])
    df_mrs["norm_reach"] = _minmax(df_mrs["cascade_reach"])
    df_mrs["norm_stab"] = _minmax(df_mrs["stability_mean"])
    df_mrs["norm_rewire"] = _minmax(df_mrs["rewiring_magnitude"])
    
    # Weighted composite score
    df_mrs["mrs_score"] = (
        0.35 * df_mrs["norm_out"] +
        0.30 * df_mrs["norm_reach"] +
        0.20 * df_mrs["norm_stab"] +
        0.15 * df_mrs["norm_rewire"]
    )
    
    df_mrs = df_mrs.sort_values(by="mrs_score", ascending=False).reset_index(drop=True)
    df_mrs["rank"] = df_mrs.index + 1
    
    logger.info(f"Top 5 Master Regulators: {list(df_mrs['gene'].head(5))}")
    return df_mrs

def export_network_graphml(
    stable_edges: pd.DataFrame,
    mrs_df: pd.DataFrame,
    gene_types: dict,
    output_path: str
) -> None:
    """Exports annotated network to GraphML format for Cytoscape and interactive dashboards."""
    G = build_networkx_graph(stable_edges, directed_only=False)
    mrs_dict = dict(zip(mrs_df["gene"], mrs_df["mrs_score"])) if mrs_df is not None else {}
    
    for n in G.nodes:
        G.nodes[n]["gene_type"] = gene_types.get(n, "effector")
        G.nodes[n]["mrs_score"] = float(mrs_dict.get(n, 0.0))
        G.nodes[n]["is_lncrna"] = int(gene_types.get(n, "") == "lncrna")
        
    nx.write_graphml(G, output_path)
    logger.info(f"Exported causal graph to {output_path} ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges).")
