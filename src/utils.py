"""
src/utils.py
Common utilities: logging, configuration loading, seed setting, and streaming downloads.
"""

import os
import sys
import yaml
import logging
import random
import urllib.request
import shutil
import numpy as np

def setup_logger(name: str = "aml_causal_atlas", log_file: str = None) -> logging.Logger:
    """Configures structured logger with console and optional file output."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    if not logger.handlers:
        formatter = logging.Formatter(
            fmt="[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        if log_file:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            
    return logger

def load_config(config_path: str = "config/config.yaml") -> dict:
    """Loads YAML configuration file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found at: {config_path}")
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return config

def set_seed(seed: int = 42) -> None:
    """Ensures complete reproducibility across NumPy, random, and Torch."""
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
    except ImportError:
        pass

def download_file_stream(url: str, dest_path: str, chunk_size: int = 8 * 1024 * 1024, expected_size: int = None) -> str:
    """
    Streams a large remote file directly to destination with resuming support.
    Ideal for streaming 6.4 GB h5ad directly to /tmp.
    """
    logger = setup_logger()
    os.makedirs(os.path.dirname(os.path.abspath(dest_path)), exist_ok=True)
    
    existing_size = os.path.getsize(dest_path) if os.path.exists(dest_path) else 0
    if expected_size and existing_size == expected_size:
        logger.info(f"File already completely downloaded at {dest_path} ({existing_size:,} bytes). Skipping.")
        return dest_path

    logger.info(f"Streaming from {url} -> {dest_path} (Existing: {existing_size:,} bytes)...")
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AML-Atlas/1.0"}
    if existing_size > 0:
        headers["Range"] = f"bytes={existing_size}-"
        
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp, open(dest_path, "ab" if existing_size > 0 else "wb") as f:
            total = int(resp.headers.get("Content-Length", 0)) + existing_size
            downloaded = existing_size
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if total > 0 and downloaded % (64 * 1024 * 1024) < chunk_size:
                    logger.info(f"Progress: {downloaded / (1024**3):.2f} GB / {total / (1024**3):.2f} GB ({downloaded / total * 100:.1f}%)")
        logger.info(f"Download complete: {dest_path} ({os.path.getsize(dest_path):,} bytes).")
    except Exception as e:
        logger.error(f"Error downloading {url}: {e}")
        raise e
        
    return dest_path
