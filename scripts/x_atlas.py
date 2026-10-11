"""Generic cell-type test for any CELLxGENE h5ad atlas (PREREG v1.11): aggregate counts per cell type, specificity, MAGMA, panel calibration.

    python scripts/x_atlas.py --build tran --files tran_amyg:AmyG tran_nac:NAc --label-col <obs column>
    python scripts/x_atlas.py --run tran
"""
import argparse
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm, rankdata
from statsmodels.stats.multitest import multipletests
import subprocess

import calibrate
import x_specific_component as x

ATL = Path("D:/AAN_data/atlases")
DERIVED = Path("D:/AAN_data/derived")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")
GENELOC = "D:/AAN_data/tools/ref/NCBI38.gene.loc"
HGNC = "D:/AAN_data/translational/hgnc_complete_set.txt"
RESULTS = Path("D:/AAN/results")
MIN_CELLS = 50


def rank_normal(v):
    r = rankdata(v, method="average")
    return norm.ppf((r - 0.5) / len(v))


def aggregate(path, region, label_col):
    a = ad.read_h5ad(path, backed="r")
    X = a.raw.X if a.raw is not None else a.X
    var = a.raw.var if a.raw is not None else a.var
    sym = var["feature_name"].astype(str).to_numpy() if "feature_name" in var.columns else var.index.astype(str).to_numpy()  # Ensembl ids if no names
    labels = a.obs[label_col].astype(str).to_numpy()
    cats = sorted({l for l in labels if (labels == l).sum() >= MIN_CELLS and not l.startswith("drop")})  # drop.* = doublets / low-quality clusters
    sums = {c: np.zeros(X.shape[1]) for c in cats}
    idx = {c: np.flatnonzero(labels == c) for c in cats}
    for c, ii in idx.items():
        for s in range(0, len(ii), 5000):
            blk = X[ii[s:s + 5000]]
            blk = blk.toarray() if sp.issparse(blk) else np.asarray(blk)
            sums[c] += blk.sum(0)
    out = pd.DataFrame({f"{region}__{c}": sums[c] for c in cats}, index=sym)
    return out.groupby(level=0).sum(), {c: len(idx[c]) for c in cats}


def build(name, files, label_col):
    parts = []
    for spec in files:
        f, region = spec.split(":")
        m, n = aggregate(ATL / f"{f}.h5ad", region, label_col)
        print(f, region, m.shape, "groups", m.shape[1])
        parts.append(m)
    M = pd.concat(parts, axis=1).fillna(0)
    cpm = M.div(M.sum(0), axis=1) * 1e6
    if str(cpm.index[0]).startswith("ENSG"):
        h = pd.read_csv(HGNC, sep="\t", low_memory=False, usecols=["ensembl_gene_id", "entrez_id"]).dropna().drop_duplicates("ensembl_gene_id")
        h["entrez_id"] = h["entrez_id"].astype(int)
        cpm = cpm.join(h.set_index("ensembl_gene_id").rename(columns={"entrez_id": "entrez"}), how="inner")
    else:
        loc = pd.read_csv(GENELOC, sep="\t", header=None, usecols=[0, 5], names=["entrez", "symbol"]).drop_duplicates("symbol")
        cpm = cpm.join(loc.set_index("symbol"), how="inner").dropna(subset=["entrez"])
    cpm["entrez"] = cpm["entrez"].astype(int)
    cpm = cpm.groupby("entrez").mean()
    cpm = cpm[cpm.mean(axis=1) >= 1]
    groups = list(cpm.columns)
    ids = {g: "G%03d_%s" % (i, re.sub(r"[^A-Za-z0-9]+", "", g)[:14]) for i, g in enumerate(groups)}
    spec = cpm.div(cpm.sum(axis=1), axis=0).apply(rank_normal).rename(columns=ids)
    spec["avg_all"] = np.log10(cpm.mean(axis=1) + 1e-3)
    spec.index.name = "GENE"
    spec.reset_index().to_csv(DERIVED / f"atlas_{name}_covar.txt", sep="\t", index=False, float_format="%.5f")
    pd.Series({v: k for k, v in ids.items()}).to_csv(DERIVED / f"atlas_{name}_names.tsv", sep="\t", header=False)
    print("genes", len(cpm), "groups", len(groups))


def magma(name, target):
    w = DERIVED / "magma" / target
    subprocess.run([str(MAGMA), "--gene-results", str(w / "genes.genes.raw"), "--gene-covar", str(DERIVED / f"atlas_{name}_covar.txt"),
                    "--model", "condition-hide=avg_all", "direction=pos", "--out", str(w / f"atlas_{name}")], check=True, capture_output=True)
    f = w / f"atlas_{name}.gsa.out"
    f = f if f.exists() else w / f"atlas_{name}.gsa.out.txt"
    d = pd.read_csv(f, comment="#", sep=r"\s+")
    d = d[d["TYPE"] == "COVAR"].set_index("VARIABLE").drop(index=["avg_all"], errors="ignore")
    d["z"] = calibrate.z_from_p(d["P"].to_numpy())
    return d


def run(name, targets=("decodeme_gwas_1", "decodeme_gwas_2")):
    panel = x.panel_traits()
    names = pd.read_csv(DERIVED / f"atlas_{name}_names.tsv", sep="\t", header=None, index_col=0)[1]
    with ThreadPoolExecutor(4) as pool:
        res = dict(zip(panel, pool.map(lambda t: magma(name, t), panel)))
    N = np.vstack([calibrate.robust_standardise(res[t]["z"].to_numpy()) for t in panel])
    out = {t: magma(name, t) for t in targets}
    g1 = out[targets[0]]
    g1["s"] = (calibrate.robust_standardise(g1["z"].to_numpy()) - N.mean(0)) / N.std(0)
    g1["fdr"] = multipletests(g1["P"], method="fdr_bh")[1]
    for t in targets[1:]:
        g1["p_" + t] = out[t]["P"]
    g1["group"] = [names[i] for i in g1.index]
    cols = ["group", "NGENES", "BETA", "SE", "P", "z", "fdr", "s"] + [c for c in g1.columns if c.startswith("p_")]
    g1[cols].sort_values("P").to_csv(RESULTS / f"x_ATLAS_{name}.tsv", sep="\t", float_format="%.5g")
    print(g1.sort_values("P").head(15)[["group", "z", "P", "fdr", "s"] + [c for c in g1.columns if c.startswith("p_")]].to_string())
    print("groups", len(g1), "| panel traits", len(panel), "| FDR<0.05:", int((g1["fdr"] < 0.05).sum()))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--build")
    ap.add_argument("--files", nargs="+")
    ap.add_argument("--label-col")
    ap.add_argument("--run")
    ap.add_argument("--targets", nargs="+", default=["decodeme_gwas_1", "decodeme_gwas_2"])
    a = ap.parse_args()
    if a.build:
        build(a.build, a.files, a.label_col)
    if a.run:
        run(a.run, tuple(a.targets))
