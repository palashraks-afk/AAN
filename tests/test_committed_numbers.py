"""Checks the headline numbers against the committed result files (no large data needed): `python -m pytest tests -q`."""
from pathlib import Path

import numpy as np
import pandas as pd

RES = Path(__file__).resolve().parent.parent / "results"
SIX = {"Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"}


def clusters():
    return pd.read_csv(RES / "clusters_decodeme_gwas_1.tsv", sep="\t")


def test_eleven_clusters_pass_fdr():
    assert int((clusters()["fdr_A"] < 0.05).sum()) == 11


def test_six_specific_clusters():
    c = clusters()
    got = set(c.loc[(c["fdr_A"] < 0.05) & (c["s"] >= 2) & (c["p_B"] < 0.05), "cluster_name"])
    assert got == SIX


def test_specificity_score_formula():
    c = clusters()
    s = (c["z_std"] - c["panel_mean"]) / c["panel_sd"]
    assert np.allclose(s, c["s"], atol=1e-3)


def test_independent_check_files_exist():
    for f in ("x_F1_families.tsv", "x_R1p_mecfs.tsv", "x_GTEx_S54.tsv", "x_ATLAS_tran.tsv", "x_ATLAS_cao.tsv", "x_M4_related_conditions.tsv"):
        assert (RES / f).exists(), f


def test_pain_overlap_number():
    m = pd.read_csv(RES / "x_M4_related_conditions.tsv", sep="\t").set_index("condition")
    assert m.loc["finngen_pain", "mean_z_six"] > m.loc["finngen_pain", "null_95th"]
