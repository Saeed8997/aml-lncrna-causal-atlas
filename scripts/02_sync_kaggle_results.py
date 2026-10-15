"""
scripts/02_sync_kaggle_results.py
Polls Kaggle kernel status, retrieves lightweight outputs, and validates local results.
"""

import os
import sys
import time
import pandas as pd
from kaggle.api.kaggle_api_extended import KaggleApi

def sync_results(kernel_slug: str = "metanoiam/aml-lncrna-mrna-multi-resolution-causal-atlas", target_dir: str = "results"):
    print(f"[1/3] Connecting to Kaggle API to monitor '{kernel_slug}'...")
    api = KaggleApi()
    api.authenticate()
    
    os.makedirs(target_dir, exist_ok=True)
    
    def get_status_info(status_obj):
        if hasattr(status_obj, "to_dict"):
            d = status_obj.to_dict()
        elif isinstance(status_obj, dict):
            d = status_obj
        else:
            d = {"status": str(getattr(status_obj, "status", "unknown")), "failureMessage": getattr(status_obj, "failure_message", None)}
        raw_st = str(d.get("status", "unknown")).lower()
        # strip enum prefix if present e.g. kernelworkerstatus.running -> running
        if "." in raw_st:
            raw_st = raw_st.split(".")[-1]
        fail_msg = d.get("failureMessage") or d.get("failure_message")
        return raw_st, fail_msg

    print("[2/3] Checking kernel status...", flush=True)
    status_data = api.kernels_status(kernel_slug)
    status, failure_msg = get_status_info(status_data)
    print(f"      Current status: {status}", flush=True)
    
    seen_log_count = [0]
    def show_logs():
        try:
            raw_logs = api.kernels_logs(kernel_slug)
            if not raw_logs:
                return
            if isinstance(raw_logs, str):
                import json
                log_items = json.loads(raw_logs)
            else:
                log_items = raw_logs
            if isinstance(log_items, list):
                new_items = log_items[seen_log_count[0]:]
                for it in new_items:
                    line = it.get("data", "").strip()
                    if line and not line.startswith("\r"):
                        print(f"      [Remote] {line}", flush=True)
                seen_log_count[0] = len(log_items)
        except Exception:
            pass

    show_logs()
    while status not in ["complete", "error", "cancel_acknowledged", "canceled", "cancelled"]:
        print(f"      Kernel is still {status}... Waiting 30s before next check.", flush=True)
        time.sleep(30)
        show_logs()
        status_data = api.kernels_status(kernel_slug)
        status, failure_msg = get_status_info(status_data)
        print(f"      Current status: {status}", flush=True)
    show_logs()
        
    if status != "complete":
        print(f"Error: Kernel ended with status '{status}'. Check Kaggle kernel logs for errors.", flush=True)
        if failure_msg:
            print(f"Failure message: {failure_msg}", flush=True)
        return False
        
    print(f"[3/3] Downloading completed results to '{target_dir}'...")
    api.kernels_output(kernel_slug, path=target_dir)
    print("Download completed! Inspecting downloaded artifacts:")
    
    found_files = []
    for root, _, files in os.walk(target_dir):
        for f in files:
            rel = os.path.relpath(os.path.join(root, f), target_dir)
            size = os.path.getsize(os.path.join(root, f))
            found_files.append((rel, size))
            print(f"  - {rel} ({size:,} bytes)")
            
    # Quick sanity check on tables
    mrs_path = os.path.join(target_dir, "tables", "mrs_rankings.csv")
    if os.path.exists(mrs_path):
        mrs_df = pd.read_csv(mrs_path)
        print(f"\nTop 5 Master Regulators from Run:")
        print(mrs_df[["rank", "gene", "gene_type", "mrs_score", "out_degree"]].head(5).to_string(index=False))
        
    val_path = os.path.join(target_dir, "tables", "crispri_validation.csv")
    if os.path.exists(val_path):
        val_df = pd.read_csv(val_path)
        print(f"\nCRISPRi Functional Validation Summary:")
        print(val_df[["screen_phenotype", "fold_enrichment", "fisher_pvalue", "validated_hit_symbols"]].to_string(index=False))
        
    return True

if __name__ == "__main__":
    slug = sys.argv[1] if len(sys.argv) > 1 else "metanoiam/aml-lncrna-mrna-multi-resolution-causal-atlas"
    sync_results(kernel_slug=slug)
