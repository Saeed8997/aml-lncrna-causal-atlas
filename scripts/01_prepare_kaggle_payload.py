"""
scripts/01_prepare_kaggle_payload.py
Programmatically bundles src/ modules and config into a self-contained kaggle_runner/run_kaggle_pipeline.py.
"""

import os
import sys
import base64
import json
import yaml

def prepare_payload():
    print("[1/3] Reading local src/ modules and config/config.yaml...")
    src_files = [
        "utils.py",
        "data_loading.py",
        "feature_selection.py",
        "preprocessing.py",
        "causal_discovery.py",
        "multi_resolution.py",
        "intervention.py",
        "graph_analysis.py",
        "validation.py"
    ]
    
    file_payloads = {}
    for fname in src_files:
        fpath = os.path.join("src", fname)
        with open(fpath, "r", encoding="utf-8") as f:
            file_payloads[fname] = f.read()
            
    with open("config/config.yaml", "r", encoding="utf-8") as f:
        config_text = f.read()
        
    encoded_files = {
        fname: base64.b64encode(content.encode("utf-8")).decode("ascii")
        for fname, content in file_payloads.items()
    }
    encoded_config = base64.b64encode(config_text.encode("utf-8")).decode("ascii")

    orchestrator_script = f'''"""
Kaggle Self-Contained Runner for AML lncRNA-mRNA Causal Regulatory Atlas
Runs on Kaggle GPU + Multi-Core environment.
"""

import os
import sys
import subprocess
import base64

print("[Kaggle Setup] Installing required high-performance causal and single-cell libraries...")
subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", "causal-learn", "anndata", "scanpy", "pyarrow", "networkx", "pyyaml"], check=False)

# Reconstruct modular src/ and config/ directories
os.makedirs("src", exist_ok=True)
os.makedirs("config", exist_ok=True)
os.makedirs("/kaggle/working/results/tables", exist_ok=True)
os.makedirs("/kaggle/working/results/figures", exist_ok=True)
os.makedirs("/kaggle/working/results/networks", exist_ok=True)
os.makedirs("/kaggle/working/results/logs", exist_ok=True)

ENCODED_FILES = {json.dumps(encoded_files)}
ENCODED_CONFIG = "{encoded_config}"

with open("config/config.yaml", "w", encoding="utf-8") as f:
    f.write(base64.b64decode(ENCODED_CONFIG).decode("utf-8"))

for fname, b64_str in ENCODED_FILES.items():
    with open(os.path.join("src", fname), "w", encoding="utf-8") as f:
        f.write(base64.b64decode(b64_str).decode("utf-8"))

with open("src/__init__.py", "w") as f:
    f.write("# Init\\n")

sys.path.insert(0, os.getcwd())

# Execute Pipeline
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.utils import set_seed, setup_logger
from src.data_loading import ensure_atlas_downloaded, load_backed_atlas, load_crispri_screen_data, identify_crispri_target_lncrnas
from src.feature_selection import select_hybrid_feature_set
from src.preprocessing import construct_metacell_atlas
from src.multi_resolution import run_multi_resolution_causal_discovery, compute_differential_causal_rewiring
from src.intervention import run_insilico_knockout_campaign
from src.graph_analysis import compute_master_regulator_scores, export_network_graphml
from src.validation import run_crispri_enrichment_validation

logger = setup_logger("kaggle_pipeline", log_file="/kaggle/working/results/logs/pipeline_execution.log")

def main():
    logger.info("=== Starting AML lncRNA Multi-Resolution Causal Atlas Pipeline ===")
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)
        
    set_seed(config["project"].get("random_seed", 42))
    
    # 1. Download/Stream to /tmp
    h5ad_path = ensure_atlas_downloaded(config)
    crispr_df = load_crispri_screen_data(config)
    crispr_targets = identify_crispri_target_lncrnas(crispr_df)
    
    # 2. Backed Mount & Feature Slicing
    adata_backed = load_backed_atlas(h5ad_path)
    features = select_hybrid_feature_set(adata_backed, config, crispr_targets=crispr_targets)
    
    logger.info(f"Slicing {{len(features)}} candidate features from {{adata_backed.n_obs:,}} cells into memory...")
    feature_indices = [adata_backed.var_names.get_loc(f) for f in features if f in adata_backed.var_names]
    
    max_cells = min(adata_backed.n_obs, 80000)
    subsample_idx = np.random.choice(adata_backed.n_obs, size=max_cells, replace=False)
    adata_sliced = adata_backed[subsample_idx, feature_indices].to_memory()
    
    # Map symbols to Ensembl IDs for robust typing
    var_symbols = {{}}
    for idx, row in adata_sliced.var.iterrows():
        sym = row.get("feature_name", idx) if "feature_name" in adata_sliced.var.columns else idx
        var_symbols[idx] = str(sym)
        
    tfs_upper = [x.upper() for x in config["features"]["key_tfs"]]
    gene_types = {{}}
    for g in adata_sliced.var_names:
        sym = var_symbols.get(g, str(g)).upper()
        if any(sym.startswith(p) for p in ["LINC", "MIR", "CATG", "HOTAIR", "NEAT", "MALAT", "AC0", "AL0", "PVT", "DANCR", "H19"]) or "-AS" in sym:
            gene_types[g] = "lncrna"
        elif sym in tfs_upper:
            gene_types[g] = "tf"
        else:
            gene_types[g] = "effector"
            
    # 3. GPU-Accelerated Metacell Aggregation
    meta_adata = construct_metacell_atlas(adata_sliced, config)
    
    # Save compact metacell h5ad (<50MB)
    try:
        meta_adata.write("/kaggle/working/results/metacell_atlas_compact.h5ad")
        logger.info("Saved compact metacell AnnData object (<50MB).")
    except Exception as e:
        logger.warning(f"Could not write h5ad: {{e}}")

    # 4. Multi-Resolution Causal Discovery
    meta_df = pd.DataFrame(meta_adata.X, columns=meta_adata.var_names)
    multi_res = run_multi_resolution_causal_discovery(meta_adata, gene_types, config)
    global_edges = multi_res["global"]
    lsc_edges = multi_res["lsc"]
    mono_edges = multi_res["mono"]
    
    # Convert edge IDs to symbols for clarity in exports
    for df in [global_edges, lsc_edges, mono_edges]:
        df["source_symbol"] = df["source"].map(lambda x: var_symbols.get(x, x))
        df["target_symbol"] = df["target"].map(lambda x: var_symbols.get(x, x))
        
    # Save edge tables
    global_edges.to_parquet("/kaggle/working/results/tables/stable_edges.parquet", index=False)
    global_edges.to_csv("/kaggle/working/results/tables/stable_edges.csv", index=False)
    lsc_edges.to_parquet("/kaggle/working/results/tables/lsc_edges.parquet", index=False)
    lsc_edges.to_csv("/kaggle/working/results/tables/lsc_edges.csv", index=False)
    mono_edges.to_parquet("/kaggle/working/results/tables/mono_edges.parquet", index=False)
    mono_edges.to_csv("/kaggle/working/results/tables/mono_edges.csv", index=False)
    
    # 5. Differential Causal Rewiring
    rewiring_df, gene_scores = compute_differential_causal_rewiring(lsc_edges, mono_edges, list(meta_adata.var_names))
    rewiring_df["source_symbol"] = rewiring_df["source"].map(lambda x: var_symbols.get(x, x))
    rewiring_df["target_symbol"] = rewiring_df["target"].map(lambda x: var_symbols.get(x, x))
    rewiring_df.to_csv("/kaggle/working/results/tables/state_rewiring.csv", index=False)
    
    # 6. SCM Fitting & In-Silico Knockouts
    ko_df = run_insilico_knockout_campaign(global_edges, meta_df, config, var_symbols=var_symbols)
    if "target_lncrna_symbol" not in ko_df.columns and "target_lncrna" in ko_df.columns:
        ko_df["target_lncrna_symbol"] = ko_df["target_lncrna"].map(lambda x: var_symbols.get(x, x))
    ko_df.to_csv("/kaggle/working/results/tables/ko_effects.csv", index=False)
    
    # 7. Master Regulator Score (MRS)
    mrs_df = compute_master_regulator_scores(global_edges, gene_scores, gene_types)
    mrs_df["symbol"] = mrs_df["gene"].map(lambda x: var_symbols.get(x, x))
    mrs_df.to_csv("/kaggle/working/results/tables/mrs_rankings.csv", index=False)
    
    # 8. Export GraphML Networks
    export_network_graphml(global_edges, mrs_df, gene_types, "/kaggle/working/results/networks/global_dag.graphml")
    export_network_graphml(lsc_edges, mrs_df, gene_types, "/kaggle/working/results/networks/lsc_dag.graphml")
    export_network_graphml(mono_edges, mrs_df, gene_types, "/kaggle/working/results/networks/mono_dag.graphml")
    
    # 9. CRISPRi Statistical Validation
    validation_hits = config["features"]["crispri_validation_hits"]
    # Check both gene IDs and symbols in validation
    mrs_val_df = mrs_df.copy()
    mrs_val_df["gene"] = mrs_val_df["symbol"]
    val_df = run_crispri_enrichment_validation(mrs_val_df, validation_hits, top_k_threshold=25)
    val_df.to_csv("/kaggle/working/results/tables/crispri_validation.csv", index=False)
    
    # 10. Publication Figures
    logger.info("Generating publication figures...")
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    
    # Figure 1: Top Master Regulators Barplot
    fig, ax = plt.subplots(figsize=(10, 6))
    top_mrs = mrs_df[mrs_df["gene_type"] == "lncrna"].head(15)
    sns.barplot(data=top_mrs, x="mrs_score", y="symbol", palette="viridis", ax=ax)
    ax.set_title("Top 15 Causal lncRNA Master Regulators in AML (MRS)", fontsize=14, fontweight="bold")
    ax.set_xlabel("Master Regulator Score (MRS)", fontsize=12)
    plt.tight_layout()
    plt.savefig("/kaggle/working/results/figures/mrs_top_lncrnas.png", dpi=300)
    plt.close()

    # Figure 2: CRISPRi Validation Enrichment
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(data=val_df, x="screen_phenotype", y="fold_enrichment", palette="magma", ax=ax)
    ax.axhline(1.0, color="gray", linestyle="--", label="Expected Random (1.0x)")
    ax.set_title("CRISPRi Hit Enrichment Among Top Causal lncRNAs", fontsize=14, fontweight="bold")
    ax.set_ylabel("Fold Enrichment (Observed / Expected)", fontsize=12)
    ax.legend()
    plt.tight_layout()
    plt.savefig("/kaggle/working/results/figures/crispri_enrichment_barplot.png", dpi=300)
    plt.close()

    # Figure 3: In-Silico Knockout Impact on Venetoclax Axis
    ven_ko = ko_df[ko_df["axis"] == "Venetoclax_Axis"]
    if not ven_ko.empty:
        pivot_ko = ven_ko.pivot_table(index="target_lncrna_symbol", columns="gene", values="percent_change", aggfunc="mean").fillna(0)
        fig, ax = plt.subplots(figsize=(9, 7))
        sns.heatmap(pivot_ko, cmap="vlag", center=0, annot=True, fmt=".1f", ax=ax)
        ax.set_title("do(lncRNA=0) Causal Impact on Venetoclax Axis (% Change)", fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig("/kaggle/working/results/figures/venetoclax_ko_heatmap.png", dpi=300)
        plt.close()

    logger.info("=== Pipeline execution completed successfully! All artifacts generated. ===")

if __name__ == "__main__":
    main()
'''

    os.makedirs("kaggle_runner", exist_ok=True)
    runner_dest = os.path.join("kaggle_runner", "run_kaggle_pipeline.py")
    with open(runner_dest, "w", encoding="utf-8") as f:
        f.write(orchestrator_script)
        
    print(f"[2/3] Generated self-contained payload: {runner_dest} ({os.path.getsize(runner_dest):,} bytes).")
    print("[3/3] Ready to deploy to Kaggle.")

if __name__ == "__main__":
    prepare_payload()
