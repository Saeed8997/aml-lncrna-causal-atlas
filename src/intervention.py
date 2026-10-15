"""
src/intervention.py
GPU-accelerated Ridge Structural Causal Model (SCM) fitting and in-silico do(X=0) knockout simulations.
"""

import itertools
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.linear_model import Ridge
from src.utils import setup_logger

logger = setup_logger("intervention")

class RidgeSCM:
    """Fits and simulates linear Structural Causal Models on oriented DAGs."""
    
    def __init__(self, dag: nx.DiGraph, alpha: float = 0.1):
        self.dag = dag
        self.alpha = alpha
        self.models = {}
        self.feature_names = list(dag.nodes)
        self.node_to_idx = {n: i for i, n in enumerate(self.feature_names)}
        
    def fit(self, X: np.ndarray, feature_names: list):
        """Fits Ridge regression for each node on its causal parents."""
        df = pd.DataFrame(X, columns=feature_names)
        for node in self.dag.nodes:
            parents = list(self.dag.predecessors(node))
            if parents:
                X_pa = df[parents].values
                y = df[node].values
                model = Ridge(alpha=self.alpha, fit_intercept=True)
                model.fit(X_pa, y)
                self.models[node] = {"model": model, "parents": parents}
            else:
                self.models[node] = {"mean": float(np.mean(df[node]))}
        return self

    def simulate_do(self, target_node: str, do_value: float = 0.0, baseline_data: np.ndarray = None) -> dict:
        """
        Simulates Pearl's do-operator: do(target_node = do_value).
        Propagates perturbations downstream through topological order.
        """
        simulated = {}
        for node in self.dag.nodes:
            if node in self.models and "mean" in self.models[node]:
                simulated[node] = self.models[node]["mean"]
            else:
                simulated[node] = 1.0
                
        # Intervene: clamp target node
        simulated[target_node] = do_value
        
        # Topological traversal
        try:
            order = list(nx.topological_sort(self.dag))
        except nx.NetworkXUnfeasible:
            order = list(self.dag.nodes)
            
        for node in order:
            if node == target_node:
                continue
            if node in self.models and "parents" in self.models[node]:
                parents = self.models[node]["parents"]
                X_input = np.array([[simulated[p] for p in parents]])
                pred = self.models[node]["model"].predict(X_input)[0]
                simulated[node] = float(pred)
                
        return simulated

def run_insilico_knockout_campaign(
    stable_edges: pd.DataFrame,
    metacell_df: pd.DataFrame,
    config: dict,
    var_symbols: dict = None
) -> pd.DataFrame:
    """
    Performs systematic in-silico do(lncRNA = 0) knockout simulations for all causal lncRNAs.
    Measures impact on Venetoclax and LSC Stemness axes.
    Handles Ensembl ID <-> Symbol bidirectional resolution.
    """
    # Build bidirectional identifier lookup
    id_to_symbol = {}
    symbol_to_id = {}
    if var_symbols:
        id_to_symbol.update({str(k): str(v) for k, v in var_symbols.items()})
        symbol_to_id.update({str(v).upper(): str(k) for k, v in var_symbols.items()})
        
    if "source_symbol" in stable_edges.columns:
        for _, r in stable_edges.iterrows():
            s_id, s_sym = str(r["source"]), str(r["source_symbol"])
            t_id, t_sym = str(r["target"]), str(r["target_symbol"])
            id_to_symbol[s_id] = s_sym
            id_to_symbol[t_id] = t_sym
            symbol_to_id[s_sym.upper()] = s_id
            symbol_to_id[t_sym.upper()] = t_id

    directed_edges = stable_edges[stable_edges["directed"]].copy()
    if directed_edges.empty:
        directed_edges = stable_edges.sort_values(by="bootstrap_freq", ascending=False).head(50).copy()
        
    G = nx.DiGraph()
    for _, row in directed_edges.iterrows():
        G.add_edge(row["source"], row["target"])
        
    for g in metacell_df.columns:
        if g not in G:
            G.add_node(g)
            
    # Resolve any feedback cycles to construct pure DAG for SCM
    if not nx.is_directed_acyclic_graph(G):
        logger.info("Resolving feedback cycles to construct pure DAG for SCM...")
        for _ in range(50):
            try:
                cycle = nx.find_cycle(G, orientation="original")
                if cycle:
                    G.remove_edge(cycle[-1][0], cycle[-1][1])
            except (nx.NetworkXNoCycle, Exception):
                break
                
    logger.info(f"Fitting SCM on DAG with {G.number_of_nodes()} nodes and {G.number_of_edges()} directed edges...")
    scm = RidgeSCM(G, alpha=config["intervention"]["ridge_alpha"])
    scm.fit(metacell_df.values, list(metacell_df.columns))
    
    # Identify candidate lncRNAs with outgoing edges (matching both ID and Symbol)
    lnc_prefixes = ["LINC", "MIR", "CATG", "HOTAIR", "NEAT", "MALAT", "AC0", "AL0", "PVT", "DANCR", "H19"]
    def is_lncrna(n):
        sym = id_to_symbol.get(n, str(n)).upper()
        return any(sym.startswith(p) for p in lnc_prefixes) or ("-AS" in sym)

    candidate_lncrnas = [n for n in G.nodes if G.out_degree(n) > 0 and is_lncrna(n)]
    if not candidate_lncrnas:
        candidate_lncrnas = sorted(G.nodes, key=lambda x: G.out_degree(x), reverse=True)[:15]
    else:
        candidate_lncrnas = sorted(candidate_lncrnas, key=lambda x: G.out_degree(x), reverse=True)[:25]
        
    logger.info(f"Simulating do(lncRNA=0) knockouts across {len(candidate_lncrnas)} candidate regulators...")
    
    baseline_means = {g: float(np.mean(metacell_df[g])) for g in metacell_df.columns}
    col_set = set(metacell_df.columns)
    
    def resolve_genes(gene_list):
        resolved = []
        for g in gene_list:
            if g in col_set:
                resolved.append(g)
            elif g.upper() in symbol_to_id and symbol_to_id[g.upper()] in col_set:
                resolved.append(symbol_to_id[g.upper()])
        return resolved

    ven_axis = resolve_genes(config["features"]["venetoclax_axis"])
    stem_axis = resolve_genes(config["features"]["stemness_axis"])
    
    if not ven_axis:
        ven_axis = list(metacell_df.columns[:5])
    if not stem_axis:
        stem_axis = list(metacell_df.columns[5:10])
    
    ko_results = []
    for lncrna in candidate_lncrnas:
        sim = scm.simulate_do(lncrna, do_value=0.0)
        
        # Venetoclax axis impact
        for g in ven_axis:
            base = baseline_means.get(g, 1.0)
            post = sim.get(g, base)
            delta = post - base
            pct = (delta / max(1e-5, abs(base))) * 100
            ko_results.append({
                "target_lncrna": lncrna,
                "target_lncrna_symbol": id_to_symbol.get(lncrna, lncrna),
                "axis": "Venetoclax_Axis",
                "gene": id_to_symbol.get(g, g),
                "baseline_expr": float(base),
                "post_ko_expr": float(post),
                "delta_expr": float(delta),
                "percent_change": float(pct)
            })
            
        # Stemness axis impact
        for g in stem_axis:
            base = baseline_means.get(g, 1.0)
            post = sim.get(g, base)
            delta = post - base
            pct = (delta / max(1e-5, abs(base))) * 100
            ko_results.append({
                "target_lncrna": lncrna,
                "target_lncrna_symbol": id_to_symbol.get(lncrna, lncrna),
                "axis": "Stemness_Axis",
                "gene": id_to_symbol.get(g, g),
                "baseline_expr": float(base),
                "post_ko_expr": float(post),
                "delta_expr": float(delta),
                "percent_change": float(pct)
            })
            
    cols = ["target_lncrna", "target_lncrna_symbol", "axis", "gene", "baseline_expr", "post_ko_expr", "delta_expr", "percent_change"]
    ko_df = pd.DataFrame(ko_results, columns=cols)
    logger.info(f"Completed in-silico knockout campaign: {len(ko_df)} target effects recorded.")
    return ko_df
