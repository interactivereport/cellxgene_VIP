"""Regression test for the default VIP DEG path, run in the VIP environment."""

import csv
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

import VIPInterface


class WelchTest(unittest.TestCase):
    def test_explicit_groups_do_not_create_expression_subsets(self):
        adapter = SimpleNamespace(
            data=SimpleNamespace(shape=(22, 2), var=pd.DataFrame({"symbol": ["first", "second"]})),
            compute_diffexp_ttest=Mock(return_value={"positive": [[0, 2.0, .01, .02]], "negative": []}),
        )
        with TemporaryDirectory() as temp:
            data = {"DEmethod": "default", "genes": [], "cellN": 22, "comGrp": ["A", "B"],
                    "cells": {"A": dict(enumerate(range(10))), "B": dict(enumerate(range(10, 22)))},
                    "data_adapter": adapter, "var_index": "symbol", "topN": "All", "CLItmp": temp,
                    "Rlib": "", "figOpt": {"img": "png", "fontsize": 10, "dpi": 100, "vectorFriendly": "No"},
                    "logFC": 10, "sigFDR": .05, "sigFC": 1.5, "labelSize": 1, "dotSize": .5,
                    "ymin": 0, "ymax": 20, "gsea": {"enable": False}}
            with patch.object(VIPInterface, "createData", side_effect=AssertionError("unexpected expression slice")), \
                 patch.object(VIPInterface.subprocess, "run", return_value=SimpleNamespace(stdout=b"plot", stderr=b"")):
                table, image, _, _ = json.loads(VIPInterface.DEG(data))
            self.assertEqual(list(Path(temp).iterdir()), [])
        mask_a, mask_b, top_n, cutoff = adapter.compute_diffexp_ttest.call_args.args
        np.testing.assert_array_equal(mask_a, np.arange(22) < 10)
        np.testing.assert_array_equal(mask_b, np.arange(22) >= 10)
        self.assertEqual((top_n, cutoff), (1, 0.01))
        self.assertEqual(next(csv.DictReader(io.StringIO(table)))["gene"], "first")
        self.assertEqual(image, "plot")


if __name__ == "__main__":
    unittest.main()
