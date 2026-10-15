"""
tests/test_causal_pipeline.py
Fast local unit test on synthetic 50-gene data to verify the entire causal inference pipeline.
"""

import unittest
import numpy as np
import pandas as pd
from src.utils import set_seed, load_config
from src.preprocessing import construct_metacell_atlas

try:
    import anndata as ad
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False
    from src.data_loading import MockAnnData

from src.causal_discovery import build_biological_background_knowledge, run_bootstrap_pc_discovery
from src.multi_resolution import run_multi_resolution_causal_discovery, compute_differential_causal_rewiring
from src.intervention import run_insilico_knockout_campaign
from src.graph_analysis import compute_master_regulator_scores
from src.validation import run_crispri_enrichment_validation

class TestCausalPipeline(unittest.TestCase):
    
    def setUp(self):
        set_seed(42)
        self.config = load_config("config/config.yaml")
        # Override parameters for lightweight testing
        self.config["causal_discovery"]["n_bootstraps"] = 3
        self.config["causal_discovery"]["stability_threshold"] = 0.30
        self.config["metacell_preprocessing"]["target_metacell_size"] = 10
        self.config["metacell_preprocessing"]["use_gpu"] = False
        
        # Define synthetic gene panel
        self.tfs = ["RUNX1", "SPI1", "MYC", "CEBPA"]
        self.lncrnas = ["MIR17HG", "AC009299.3", "HOTAIR", "LINC001", "LINC002", "CATG001"]
        self.effectors = ["BCL2", "MCL1", "BAX", "CD34", "CD14", "LYZ"]
        self.genes = self.tfs + self.lncrnas + self.effectors
        
        self.gene_types = {}
        for g in self.tfs:
            self.gene_types[g] = "tf"
        for g in self.lncrnas:
            self.gene_types[g] = "lncrna"
        for g in self.effectors:
            self.gene_types[g] = "effector"
            
        # Create synthetic single-cell expression (120 cells x 16 genes)
        n_cells = 120
        # True linear causal generation: TF -> lncRNA -> Effector
        np_tfs = np.random.normal(2.0, 0.5, (n_cells, len(self.tfs)))
        
        W_tf_lnc = np.random.normal(0.4, 0.1, (len(self.tfs), len(self.lncrnas)))
        np_lncs = np_tfs @ W_tf_lnc + np.random.normal(1.0, 0.2, (n_cells, len(self.lncrnas)))
        
        W_lnc_eff = np.random.normal(0.5, 0.1, (len(self.lncrnas), len(self.effectors)))
        np_effs = np_lncs @ W_lnc_eff + np.random.normal(1.2, 0.2, (n_cells, len(self.effectors)))
        
        X = np.hstack([np_tfs, np_lncs, np_effs])
        
        obs = pd.DataFrame({
            "donor_id": ["Donor1"] * 60 + ["Donor2"] * 60,
            "cell_type": (["LSC_prog"] * 30 + ["Monocytic_blast"] * 30) * 2
        })
        var = pd.DataFrame(index=self.genes)
        var["feature_name"] = self.genes
        var["feature_biotype"] = [
            "protein_coding" if g in self.tfs + self.effectors else "lncRNA" 
            for g in self.genes
        ]
        
        if HAS_ANNDATA:
            self.synth_adata = ad.AnnData(X=X, obs=obs, var=var)
        else:
            self.synth_adata = MockAnnData(X=X, obs=obs, var=var)

    def test_01_metacell_construction(self):
        meta = construct_metacell_atlas(self.synth_adata, self.config)
        self.assertGreater(meta.n_obs, 2)
        self.assertEqual(meta.n_vars, len(self.genes))
        self.assertIn("macro_state", meta.obs.columns)

    def test_02_background_knowledge(self):
        bk = build_biological_background_knowledge(self.genes, self.gene_types)
        self.assertIsNotNone(bk)

    def test_03_causal_discovery_and_rewiring(self):
        meta = construct_metacell_atlas(self.synth_adata, self.config)
        multi_res = run_multi_resolution_causal_discovery(meta, self.gene_types, self.config)
        self.assertIn("global", multi_res)
        self.assertIn("lsc", multi_res)
        self.assertIn("mono", multi_res)
        
        rewiring_df, gene_scores = compute_differential_causal_rewiring(
            multi_res["lsc"], multi_res["mono"], self.genes
        )
        self.assertIsNotNone(rewiring_df)
        self.assertEqual(len(gene_scores), len(self.genes))

    def test_04_intervention_and_validation(self):
        meta = construct_metacell_atlas(self.synth_adata, self.config)
        meta_df = pd.DataFrame(meta.X, columns=meta.var_names)
        
        multi_res = run_multi_resolution_causal_discovery(meta, self.gene_types, self.config)
        global_edges = multi_res["global"]
        
        # Run in-silico knockout
        ko_df = run_insilico_knockout_campaign(global_edges, meta_df, self.config)
        self.assertFalse(ko_df.empty)
        self.assertIn("percent_change", ko_df.columns)
        
        # Master Regulator Score
        _, gene_scores = compute_differential_causal_rewiring(
            multi_res["lsc"], multi_res["mono"], self.genes
        )
        mrs_df = compute_master_regulator_scores(global_edges, gene_scores, self.gene_types)
        self.assertEqual(len(mrs_df), len(self.genes))
        self.assertIn("mrs_score", mrs_df.columns)
        
        # CRISPRi validation
        crispri_mock = {
            "proliferation": ["MIR17HG", "HOTAIR"],
            "venetoclax": ["AC009299.3"]
        }
        val_df = run_crispri_enrichment_validation(mrs_df, crispri_mock, top_k_threshold=5)
        self.assertEqual(len(val_df), 2)
        self.assertIn("fold_enrichment", val_df.columns)

if __name__ == "__main__":
    unittest.main()
