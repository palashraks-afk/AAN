"""PREREG v1.9 R1': cell-type test for ME/CFS conditioning on the principal components of the other panel traits (MAGMA).

    python scripts/x_pc_conditioning.py --validate     # gate V' on controls only
    python scripts/x_pc_conditioning.py --mecfs        # panel null, ME/CFS gwas_1 and gwas_2, calls
"""
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import calibrate
import x_specific_component as x

DERIVED = Path("D:/AAN_data/derived")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")
RESULTS = Path("D:/AAN/results")


def run_target(target, others, cov):
    genes = x.common_genes(cov, others + [target])
    Z = np.column_stack([x.load_genes(o).loc[genes, "ZSTAT"].to_numpy() for o in others])
    pcs, k, expl = x.kaiser_pcs(Z)
    df = cov.loc[genes].copy()
    for i in range(k):
        df[f"PC{i + 1}"] = pcs[:, i]
    workdir = DERIVED / "magma" / target
    covar = workdir / f"pc_covar_{len(others)}.txt"
    df.reset_index().to_csv(covar, sep="\t", index=False, float_format="%.5f")
    out = workdir / ("PCall" if target.startswith("decodeme") else "PCloo")
    cond = ["avg_all", "avg_neuron"] + [f"PC{i + 1}" for i in range(k)]
    subprocess.run([str(MAGMA), "--gene-results", str(workdir / "genes.genes.raw"), "--gene-covar", str(covar),
                    "--model", "condition-hide=" + ",".join(cond), "direction=pos", "--out", str(out)],
                   check=True, capture_output=True)
    f = Path(str(out) + ".gsa.out")
    f = f if f.exists() else Path(str(out) + ".gsa.out.txt")
    df = pd.read_csv(f, comment="#", sep=r"\s+")
    df = df[df["VARIABLE"].str.match(r"^c\d+$")].copy()
    df["cluster_id"] = df["VARIABLE"].str[1:].astype(int)
    df = df.set_index("cluster_id").sort_index()
    df["z"] = calibrate.z_from_p(df["P"].to_numpy())
    df["fdr"] = multipletests(df["P"], method="fdr_bh")[1]
    df["k"] = k
    return df


def validate():
    cov, _ = x.load_covar()
    panel = x.panel_traits()
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    rows = []
    a = run_target("alzheimer", [p for p in panel if p != "alzheimer"], cov).join(ann)
    v2 = bool(((a["supercluster"] == "Microglia") & (a["fdr"] < 0.05)).any())
    print("V2' alzheimer keeps a microglia cluster at FDR<0.05:", v2, "| k =", int(a["k"].iloc[0]))
    h = run_target("height", [p for p in panel if p != "height"], cov).join(ann)
    top3 = h.sort_values("z", ascending=False).head(3)
    v3 = bool(top3["supercluster"].isin(calibrate.NON_NEURONAL).all())
    print("V3' height top-3:", top3["supercluster"].tolist(), "non-neuronal:", v3)
    pd.DataFrame([{"check": "V2_alzheimer_microglia", "pass": v2}, {"check": "V3_height_top3_nonneuronal", "pass": v3}]) \
        .to_csv(RESULTS / "x_R1p_validation.tsv", sep="\t", index=False)
    print("GATE V':", "PASS" if v2 and v3 else "FAIL")
    return v2 and v3


def mecfs():
    cov, _ = x.load_covar()
    panel = x.panel_traits()
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    with ThreadPoolExecutor(4) as pool:
        res = list(pool.map(lambda t: run_target(t, [p for p in panel if p != t], cov), panel))
    N = np.vstack([calibrate.robust_standardise(r["z"].to_numpy()) for r in res])
    out = {}
    for name in ("decodeme_gwas_1", "decodeme_gwas_2"):
        r = run_target(name, panel, cov)
        r["s"] = (calibrate.robust_standardise(r["z"].to_numpy()) - N.mean(0)) / N.std(0)
        out[name] = r
    g1, g2 = out["decodeme_gwas_1"], out["decodeme_gwas_2"]
    t = g1[["BETA", "SE", "P", "z", "fdr", "s", "k"]].join(g2[["P", "z", "s"]], rsuffix="_gwas2").join(ann[["cluster_name", "supercluster"]])
    t["called"] = (t["fdr"] < 0.05) & (t["P_gwas2"] < 0.05) & (t["s"] >= 2.0)
    t.to_csv(RESULTS / "x_R1p_mecfs.tsv", sep="\t", float_format="%.5g")
    print("panel traits:", len(panel), "| k =", int(g1["k"].iloc[0]))
    print("gwas_1 clusters FDR<0.05 after conditioning on panel PCs:", int((g1["fdr"] < 0.05).sum()))
    print("called ME/CFS-specific:", int(t["called"].sum()))
    print(t.sort_values("P").head(12)[["cluster_name", "supercluster", "z", "fdr", "s", "P_gwas2", "called"]].to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--mecfs", action="store_true")
    a = ap.parse_args()
    ok = True
    if a.validate:
        ok = validate()
    if a.mecfs:
        if not ok:
            raise SystemExit("gate V' failed: R1' is not run")
        mecfs()
