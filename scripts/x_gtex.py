"""PREREG v1.10: GTEx bulk-tissue check of the amygdala finding (H1), peripheral nerve (H2) and pituitary (H3).

    python scripts/x_gtex.py --build
    python scripts/x_gtex.py --run
"""
import argparse
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
from statsmodels.stats.multitest import multipletests

import calibrate
import x_specific_component as x

DERIVED = Path("D:/AAN_data/derived")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")
GTEX = Path("D:/AAN_data/translational/gtex_median_tpm.gct.gz")
GENELOC = "D:/AAN_data/tools/ref/NCBI38.gene.loc"
RESULTS = Path("D:/AAN/results")


def clean(t):
    return re.sub(r"[^A-Za-z0-9]+", "_", t).strip("_")


def short(t, tissues):
    return "T%02d_%s" % (tissues.index(t), clean(t)[:12])  # MAGMA cuts long variable names


def rank_normal(v):
    r = rankdata(v, method="average")
    return norm.ppf((r - 0.5) / len(v))


def build():
    g = pd.read_csv(GTEX, sep="\t", skiprows=2)
    loc = pd.read_csv(GENELOC, sep="\t", header=None, usecols=[0, 5], names=["entrez", "symbol"])
    g = g.merge(loc.drop_duplicates("symbol"), left_on="Description", right_on="symbol")
    tissues = [c for c in g.columns if c not in ("Name", "Description", "entrez", "symbol")]
    g = g.drop_duplicates("entrez").set_index("entrez")[tissues]
    g = g[g.mean(axis=1) >= 1]
    brain = [t for t in tissues if t.startswith("Brain")]
    print("tissues", len(tissues), "brain", len(brain), "genes", len(g))
    avg54 = np.log10(g.mean(axis=1) + 1e-3)
    avgb = np.log10(g[brain].mean(axis=1) + 1e-3)
    print("genes with brain expression:", int((g[brain].sum(axis=1) > 0).sum()))
    s54 = g.div(g.sum(axis=1), axis=0)
    gb = g[g[brain].sum(axis=1) > 0]  # within-brain fraction is undefined for genes with no brain expression
    sb = gb[brain].div(gb[brain].sum(axis=1), axis=0)
    ids = {t: short(t, tissues) for t in tissues}
    S54 = s54.apply(rank_normal).rename(columns=ids)
    SB = sb.apply(rank_normal).rename(columns=ids)
    for df, name in ((S54, "S54"), (SB, "SB")):
        out = df.copy()
        out["avg54"] = avg54.loc[out.index]
        out["avgbrain"] = avgb.loc[out.index]
        out.index.name = "GENE"
        out.reset_index().to_csv(DERIVED / f"gtex_{name}_covar.txt", sep="\t", index=False, float_format="%.5f")
    pd.Series({ids[t]: t for t in tissues}).to_csv(DERIVED / "gtex_tissue_names.tsv", sep="\t", header=False)


def magma(target, which):
    workdir = DERIVED / "magma" / target
    cond = "avg54" if which == "S54" else "avg54,avgbrain"
    subprocess.run([str(MAGMA), "--gene-results", str(workdir / "genes.genes.raw"), "--gene-covar",
                    str(DERIVED / f"gtex_{which}_covar.txt"), "--model", "condition-hide=" + cond, "direction=pos",
                    "--out", str(workdir / f"gtex_{which}")], check=True, capture_output=True)
    f = workdir / f"gtex_{which}.gsa.out"
    f = f if f.exists() else workdir / f"gtex_{which}.gsa.out.txt"
    d = pd.read_csv(f, comment="#", sep=r"\s+")
    d = d[d["TYPE"] == "COVAR"].set_index("VARIABLE")
    d = d.drop(index=[i for i in ("avg54", "avgbrain") if i in d.index])
    d["z"] = calibrate.z_from_p(d["P"].to_numpy())
    return d


def run():
    panel = x.panel_traits()
    names = pd.read_csv(DERIVED / "gtex_tissue_names.tsv", sep="\t", header=None, index_col=0)[1]
    for which in ("S54", "SB"):
        with ThreadPoolExecutor(4) as pool:
            res = dict(zip(panel, pool.map(lambda t: magma(t, which), panel)))
        N = np.vstack([calibrate.robust_standardise(res[t]["z"].to_numpy()) for t in panel])
        g1, g2 = magma("decodeme_gwas_1", which), magma("decodeme_gwas_2", which)
        g1["s"] = (calibrate.robust_standardise(g1["z"].to_numpy()) - N.mean(0)) / N.std(0)
        g1["fdr"] = multipletests(g1["P"], method="fdr_bh")[1]
        g1["p_gwas2"] = g2["P"]
        g1["tissue"] = [names[i] for i in g1.index]
        g1[["tissue", "NGENES", "BETA", "SE", "P", "z", "fdr", "p_gwas2", "s"]].sort_values("P") \
            .to_csv(RESULTS / f"x_GTEx_{which}.tsv", sep="\t", float_format="%.5g")
        print(which, "panel traits", len(panel))
        print(g1.sort_values("P").head(10)[["tissue", "z", "P", "fdr", "p_gwas2", "s"]].to_string())
        if which == "S54":
            for t in ("Nerve - Tibial", "Pituitary"):
                print(t, g1[g1["tissue"] == t][["P", "p_gwas2", "s"]].to_dict("records"))
        else:
            print("Amygdala", g1[g1["tissue"] == "Brain - Amygdala"][["P", "p_gwas2", "s", "fdr"]].to_dict("records"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--run", action="store_true")
    a = ap.parse_args()
    if a.build:
        build()
    if a.run:
        run()
