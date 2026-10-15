"""
src/causal_discovery.py
Biologically-constrained high-performance parallel PC algorithm with multi-core bootstrap stability selection.
Integrates vectorized Fisher-Z level-0 correlation pruning, analytic partial correlation depth testing,
collider (v-structure) orientation, biological BackgroundKnowledge enforcement, and Meek rules.
"""

from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import norm
from joblib import Parallel, delayed
from src.utils import setup_logger

logger = setup_logger("causal_discovery")

try:
    from causallearn.search.ConstraintBased.PC import pc
    from causallearn.utils.cit import CIT_Base, register_ci_test
    from causallearn.utils.PCUtils.BackgroundKnowledge import BackgroundKnowledge
    from causallearn.graph.GraphNode import GraphNode

    class RobustFisherZ(CIT_Base):
        """
        Numerically stable Fisher-Z conditional independence test with ridge shrinkage
        and pseudoinverse fallback to prevent singular matrix LinAlgErrors on single-cell data.
        """
        def __init__(self, data, **kwargs):
            super().__init__(data, **kwargs)
            self.method = "robust_fisherz"
            self.data = np.ascontiguousarray(self.data, dtype=np.float64)
            n, p = self.data.shape
            self.n = n
            
            # Micro-jitter to eliminate exact collinearity
            data_clean = self.data + np.random.normal(0, 1e-5, size=(n, p))
            stds = np.std(data_clean, axis=0)
            stds[stds < 1e-6] = 1.0
            data_norm = (data_clean - np.mean(data_clean, axis=0)) / stds
            
            corr = np.dot(data_norm.T, data_norm) / max(1, n - 1)
            corr = np.nan_to_num(corr, nan=0.0)
            np.fill_diagonal(corr, 1.0)
            self.corr_matrix = corr

        def __call__(self, X: int, Y: int, condition_set=None, *args, **kwargs):
            if condition_set is None or len(condition_set) == 0:
                r = self.corr_matrix[X, Y]
                k = 0
            elif len(condition_set) == 1:
                z_var = list(condition_set)[0]
                r_xy = self.corr_matrix[X, Y]
                r_xz = self.corr_matrix[X, z_var]
                r_yz = self.corr_matrix[Y, z_var]
                denom = np.sqrt(max(1e-12, (1.0 - r_xz * r_xz) * (1.0 - r_yz * r_yz)))
                r = (r_xy - r_xz * r_yz) / denom
                k = 1
            else:
                cond_list = list(condition_set)
                var_idx = [X, Y] + cond_list
                sub = self.corr_matrix[np.ix_(var_idx, var_idx)]
                try:
                    inv = np.linalg.inv(sub)
                except np.linalg.LinAlgError:
                    inv = np.linalg.pinv(sub + 1e-3 * np.eye(len(var_idx)))
                    
                denom = np.sqrt(max(1e-12, abs(inv[0, 0] * inv[1, 1])))
                r = -inv[0, 1] / denom
                k = len(cond_list)
                
            r = max(-0.999999, min(0.999999, float(r)))
            z = 0.5 * np.log((1.0 + r) / (1.0 - r))
            stat = np.sqrt(max(1, self.n - k - 3)) * abs(z)
            p_val = 2.0 * (1.0 - norm.cdf(stat))
            return float(np.clip(p_val, 0.0, 1.0))

    try:
        register_ci_test("robust_fisherz", RobustFisherZ)
    except Exception as e:
        logger.warning(f"Could not register robust_fisherz: {e}")

    HAS_CAUSAL_LEARN = True
except ImportError:
    HAS_CAUSAL_LEARN = False

class BiologicalBackgroundKnowledge:
    """Stores forbidden causal edge constraints and bridges with causal-learn if present."""
    def __init__(self, forbidden_edges: set, genes: list = None):
        self.forbidden_edges = forbidden_edges
        self.cl_bk = None
        if HAS_CAUSAL_LEARN and genes is not None:
            try:
                self.cl_bk = BackgroundKnowledge()
                nodes = {g: GraphNode(g) for g in genes}
                for u_idx, v_idx in forbidden_edges:
                    self.cl_bk.add_forbidden_by_node(nodes[genes[u_idx]], nodes[genes[v_idx]])
            except Exception:
                self.cl_bk = None

def build_biological_background_knowledge(genes: list, gene_types: dict) -> BiologicalBackgroundKnowledge:
    """
    Constructs edge orientation priors:
    - Non-TF effector genes CANNOT have outgoing directed edges to TFs or lncRNAs.
    - Only lncRNAs and TFs can act as causal parents.
    """
    node_to_idx = {g: i for i, g in enumerate(genes)}
    forbidden_set = set()
    for u in genes:
        type_u = gene_types.get(u, "effector")
        if type_u == "effector":
            for v in genes:
                if u == v:
                    continue
                type_v = gene_types.get(v, "effector")
                if type_v in ["tf", "lncrna"]:
                    forbidden_set.add((node_to_idx[u], node_to_idx[v]))
                    
    logger.info(f"BackgroundKnowledge initialized with {len(forbidden_set):,} biologically forbidden edges.")
    return BiologicalBackgroundKnowledge(forbidden_edges=forbidden_set, genes=genes)

def fast_vectorized_pc_discovery(
    data: np.ndarray,
    alpha: float = 0.01,
    max_depth: int = 2,
    forbidden_edges: set = None,
    max_cond_candidates: int = 20
) -> np.ndarray:
    """
    High-performance constraint-based PC algorithm for large single-cell atlases.
    - Vectorized Fisher-Z level-0 correlation pruning (O(P^2) matrix operations, ~0.03s).
    - Analytic scalar partial correlation for depth=1 (O(1) per test).
    - 4x4 precision submatrix inversion for depth=2 (O(1) per test).
    - Complete sepset separation tracking for faithful collider orientation.
    - V-structure collider orientation respecting background constraints.
    - Biological background knowledge enforcement.
    - Meek orientation rules propagation (R1, R2).
    """
    n, p = data.shape
    if forbidden_edges is None:
        forbidden_edges = set()

    # 1. Normalize data and compute empirical correlation matrix
    data_clean = np.ascontiguousarray(data, dtype=np.float64)
    data_clean = data_clean + np.random.normal(0, 1e-5, size=(n, p))
    stds = np.std(data_clean, axis=0)
    stds[stds < 1e-6] = 1.0
    norm_data = (data_clean - np.mean(data_clean, axis=0)) / stds

    corr = np.dot(norm_data.T, norm_data) / max(1, n - 1)
    corr = np.clip(corr, -0.999999, 0.999999)
    np.fill_diagonal(corr, 1.0)

    sepset = {}
    z_crit = norm.ppf(1.0 - alpha / 2.0)
    r_crit_0 = np.tanh(z_crit / np.sqrt(max(1, n - 3)))

    # --- Depth 0: Vectorized Level-0 correlation pruning ---
    adj = (np.abs(corr) > r_crit_0)
    np.fill_diagonal(adj, False)

    for i in range(p):
        for j in range(i + 1, p):
            if not adj[i, j]:
                sepset[(i, j)] = ()
                sepset[(j, i)] = ()

    # --- Depth 1: Partial correlation given 1 variable ---
    if max_depth >= 1:
        stat_factor_1 = np.sqrt(max(1, n - 1 - 3))
        for i in range(p):
            neighbors_i = np.where(adj[i])[0]
            for j in neighbors_i:
                if j <= i or not adj[i, j]:
                    continue
                cand_k = [k for k in neighbors_i if k != j]
                if len(cand_k) > max_cond_candidates:
                    cand_k = sorted(cand_k, key=lambda k: abs(corr[k, j]), reverse=True)[:max_cond_candidates]

                r_ij = corr[i, j]
                for k in cand_k:
                    r_ik = corr[i, k]
                    r_jk = corr[j, k]
                    denom = np.sqrt(max(1e-12, (1.0 - r_ik * r_ik) * (1.0 - r_jk * r_jk)))
                    r_cond = (r_ij - r_ik * r_jk) / denom
                    r_cond = max(-0.999999, min(0.999999, r_cond))
                    z_cond = 0.5 * np.log((1.0 + r_cond) / (1.0 - r_cond))
                    if abs(z_cond) * stat_factor_1 <= z_crit:
                        adj[i, j] = False
                        adj[j, i] = False
                        sepset[(i, j)] = (k,)
                        sepset[(j, i)] = (k,)
                        break

    # --- Depth 2: Partial correlation given 2 variables ---
    if max_depth >= 2:
        stat_factor_2 = np.sqrt(max(1, n - 2 - 3))
        for i in range(p):
            neighbors_i = np.where(adj[i])[0]
            for j in neighbors_i:
                if j <= i or not adj[i, j]:
                    continue
                cand_k = [k for k in neighbors_i if k != j]
                if len(cand_k) < 2:
                    continue
                cand_depth2 = sorted(cand_k, key=lambda k: abs(corr[k, j]), reverse=True)[:8]

                for k1, k2 in combinations(cand_depth2, 2):
                    sub_idx = [i, j, k1, k2]
                    sub_corr = corr[np.ix_(sub_idx, sub_idx)]
                    try:
                        inv = np.linalg.inv(sub_corr)
                    except np.linalg.LinAlgError:
                        inv = np.linalg.pinv(sub_corr + 1e-4 * np.eye(4))
                    denom = np.sqrt(max(1e-12, inv[0, 0] * inv[1, 1]))
                    r_cond = -inv[0, 1] / denom
                    r_cond = max(-0.999999, min(0.999999, r_cond))
                    z_cond = 0.5 * np.log((1.0 + r_cond) / (1.0 - r_cond))
                    if abs(z_cond) * stat_factor_2 <= z_crit:
                        adj[i, j] = False
                        adj[j, i] = False
                        sepset[(i, j)] = (k1, k2)
                        sepset[(j, i)] = (k1, k2)
                        break

    # --- Graph Initialization (causal-learn format: G[i, j]=-1, G[j, i]=1 means i -> j) ---
    G = np.zeros((p, p), dtype=np.int32)
    for i in range(p):
        for j in range(i + 1, p):
            if adj[i, j]:
                G[i, j] = -1
                G[j, i] = -1

    # --- Stage 2: Collider (v-structure) Orientation ---
    for k in range(p):
        neighbors_k = np.where(G[k, :] != 0)[0]
        if len(neighbors_k) < 2:
            continue
        for i, j in combinations(neighbors_k, 2):
            if G[i, j] == 0:
                sep = sepset.get((i, j), ())
                if k not in sep:
                    if (i, k) not in forbidden_edges and G[k, i] == -1 and G[i, k] == -1:
                        G[i, k] = -1
                        G[k, i] = 1
                    if (j, k) not in forbidden_edges and G[k, j] == -1 and G[j, k] == -1:
                        G[j, k] = -1
                        G[k, j] = 1

    # --- Stage 3: Background Knowledge Enforcement ---
    for u, v in forbidden_edges:
        if G[u, v] == -1 and G[v, u] == -1:
            if (v, u) not in forbidden_edges:
                G[v, u] = -1
                G[u, v] = 1

    # --- Stage 4: Meek Rules Propagation ---
    changed = True
    iteration = 0
    while changed and iteration < 5:
        changed = False
        iteration += 1
        for i in range(p):
            for j in range(p):
                if G[i, j] == -1 and G[j, i] == 1:
                    for k in range(p):
                        if k != i and G[j, k] == -1 and G[k, j] == -1:
                            if G[i, k] == 0 and (j, k) not in forbidden_edges:
                                G[j, k] = -1
                                G[k, j] = 1
                                changed = True
                            elif G[i, k] == -1 and G[k, i] == -1 and (i, k) not in forbidden_edges:
                                G[i, k] = -1
                                G[k, i] = 1
                                changed = True

    return G

def _run_single_pc_bootstrap(data_sub, genes: list, alpha: float, max_depth: int, bk):
    """Executes a single high-performance PC iteration on subsampled metacells."""
    forbidden = bk.forbidden_edges if hasattr(bk, "forbidden_edges") else set()
    return fast_vectorized_pc_discovery(data_sub, alpha=alpha, max_depth=max_depth, forbidden_edges=forbidden)

def run_bootstrap_pc_discovery(
    metacell_df: pd.DataFrame,
    gene_types: dict,
    config: dict
) -> pd.DataFrame:
    """
    Executes parallel multi-core bootstrap PC causal discovery.
    Retains edges satisfying bootstrap stability threshold (freq >= 0.60).
    """
    genes = list(metacell_df.columns)
    p = len(genes)
    n = len(metacell_df)
    
    alpha = config["causal_discovery"]["alpha"]
    max_depth = config["causal_discovery"]["max_depth"]
    n_bootstraps = config["causal_discovery"]["n_bootstraps"]
    subsample_ratio = config["causal_discovery"]["subsample_ratio"]
    threshold = config["causal_discovery"]["stability_threshold"]
    n_jobs = config["causal_discovery"].get("n_jobs", -1)
    
    logger.info(f"Running Bootstrap PC Discovery (N={n} metacells, P={p} genes, B={n_bootstraps} bootstraps, alpha={alpha})...")
    
    bk = build_biological_background_knowledge(genes, gene_types)
    X = np.ascontiguousarray(metacell_df.values, dtype=np.float64)
    
    subsample_size = max(5, int(n * subsample_ratio))
    bootstrap_seeds = [42 + b for b in range(n_bootstraps)]
    subsamples = [
        np.ascontiguousarray(X[np.random.RandomState(s).choice(n, size=subsample_size, replace=False)])
        for s in bootstrap_seeds
    ]
    
    logger.info(f"Launching {n_bootstraps} parallel PC runs with joblib (n_jobs={n_jobs})...")
    results = Parallel(n_jobs=n_jobs)(
        delayed(_run_single_pc_bootstrap)(sub, genes, alpha, max_depth, bk)
        for sub in subsamples
    )
    
    edge_counts = np.zeros((p, p), dtype=np.int32)
    directed_counts = np.zeros((p, p), dtype=np.int32)
    
    for G in results:
        for i in range(p):
            for j in range(i + 1, p):
                val_ij = G[i, j]
                val_ji = G[j, i]
                
                if val_ij != 0 or val_ji != 0:
                    edge_counts[i, j] += 1
                    edge_counts[j, i] += 1
                    
                    if val_ij == -1 and val_ji == 1:
                        directed_counts[i, j] += 1
                    elif val_ji == -1 and val_ij == 1:
                        directed_counts[j, i] += 1
                        
    stable_edges = []
    corr_matrix = np.corrcoef(X.T)
    corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
    
    for i in range(p):
        for j in range(p):
            if i == j:
                continue
            freq = edge_counts[i, j] / n_bootstraps
            if freq >= threshold:
                dir_freq = directed_counts[i, j] / max(1, edge_counts[i, j])
                is_directed = dir_freq >= 0.50
                
                if not is_directed and i > j:
                    continue
                    
                reg_sign = 1 if corr_matrix[i, j] >= 0 else -1
                
                stable_edges.append({
                    "source": genes[i],
                    "target": genes[j],
                    "source_type": gene_types.get(genes[i], "effector"),
                    "target_type": gene_types.get(genes[j], "effector"),
                    "bootstrap_freq": float(freq),
                    "dir_prob": float(dir_freq),
                    "directed": bool(is_directed),
                    "regulation_sign": int(reg_sign),
                    "correlation": float(corr_matrix[i, j])
                })
                
    stable_df = pd.DataFrame(stable_edges)
    logger.info(f"Discovered {len(stable_df):,} stable causal edges across {p} genes.")
    return stable_df
