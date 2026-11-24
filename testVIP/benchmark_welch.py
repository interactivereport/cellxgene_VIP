"""Measure VIP's group preparation on a deterministic synthetic expression matrix.

Run in a fresh Python process in the VIP environment from the checkout root.
Pass --module-source a saved copy of the original VIPInterface.py to compare.
The CELLxGENE test and R rendering are stubbed to isolate VIP's matrix copies.
"""

import argparse
import importlib.util
import resource
import sys
from tempfile import TemporaryDirectory
import time
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import VIPInterface


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module-source")
    args = parser.parse_args()
    vip = VIPInterface
    if args.module_source:
        spec = importlib.util.spec_from_file_location("benchmark_vip", args.module_source)
        vip = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(vip)

    rows, columns, entries = 64_000, 4096, 96
    rng = np.random.default_rng(87)
    ptr = np.arange(rows + 1, dtype=np.int32) * entries
    indices = rng.integers(columns, size=rows * entries, dtype=np.int32)
    values = rng.random(rows * entries, dtype=np.float32)
    matrix = sparse.csr_matrix((values, indices, ptr), shape=(rows, columns))
    matrix.sum_duplicates()
    obs = pd.DataFrame({"label": np.arange(rows).astype(str)}, index=[str(i) for i in range(rows)])
    var = pd.DataFrame({"symbol": [f"g{i}" for i in range(columns)]}, index=[str(i) for i in range(columns)])
    adata = ad.AnnData(matrix, obs=obs, var=var)
    adapter = SimpleNamespace(
        data=adata, get_obs_keys=lambda: ["label"],
        compute_diffexp_ttest=lambda *_: {"positive": [[0, 2.0, .01, .02]], "negative": []},
    )
    with TemporaryDirectory() as temp:
        data = {"DEmethod": "default", "genes": [], "cellN": rows, "comGrp": ["A", "B"],
                "cells": {"A": dict(enumerate(range(rows // 2))), "B": dict(enumerate(range(rows // 2, rows)))},
                "data_adapter": adapter, "obs_index": "label", "var_index": "symbol", "topN": "All",
                "CLItmp": temp, "Rlib": "", "figOpt": {"img": "png", "fontsize": 10, "dpi": 100,
                "vectorFriendly": "No", "scale": "No"}, "logFC": 10, "sigFDR": .05, "sigFC": 1.5,
                "labelSize": 1, "dotSize": .5, "ymin": 0, "ymax": 20, "gsea": {"enable": False}}
        with patch.object(vip.subprocess, "run", return_value=SimpleNamespace(stdout=b"plot", stderr=b"")):
            start = time.perf_counter()
            vip.DEG(data)
            elapsed = time.perf_counter() - start
    mib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**2 if sys.platform == "darwin" else 1024)
    print(f"{matrix.shape}, nnz={matrix.nnz}: {elapsed:.3f} s, peak RSS {mib:.1f} MiB")


if __name__ == "__main__":
    main()
