"""Cell-type results for ME/CFS and the specificity calibration against the comparison panel.

Reads the MAGMA gene-property output for every trait (models A and B), applies the rules frozen in
preregistration/PREREG_v1.md and writes tables to results/.

    python scripts/calibrate.py --target decodeme_gwas_1
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

DERIVED = Path("D:/AAN_data/derived")
RESULTS = Path("D:/AAN/results")

EXCITATORY = {
    "Upper-layer intratelencephalic", "Deep-layer intratelencephalic", "Deep-layer near-projecting",
    "Deep-layer corticothalamic and 6b",
}
NON_NEURONAL = {
    "Astrocyte", "Microglia", "Oligodendrocyte", "Oligodendrocyte precursor",
    "Committed oligodendrocyte precursor", "Vascular", "Fibroblast", "Ependymal",
    "Choroid plexus", "Bergmann glia", "Miscellaneous",
}


def load_model(name, model):
    path = DERIVED / "magma" / name / f"cluster_{model}.gsa.out"
    df = pd.read_csv(path, comment="#", sep=r"\s+")
    df = df[df["VARIABLE"].str.match(r"^c\d+$")].copy()
    df["cluster_id"] = df["VARIABLE"].str[1:].astype(int)
    return df.set_index("cluster_id").sort_index()


def z_from_p(p):
    return norm.isf(np.clip(p, 1e-300, 1 - 1e-16))


def robust_standardise(z):
    z = np.asarray(z, dtype=float)
    mad = np.median(np.abs(z - np.median(z))) * 1.4826
    return (z - np.median(z)) / mad


def traits_with_results(model="A"):
    base = DERIVED / "magma"
    return sorted(d.name for d in base.iterdir() if (d / f"cluster_{model}.gsa.out").exists())


def passes_heritability_gate(names):
    path = DERIVED / "ldsc_results.tsv"
    if not path.exists():
        return {n: None for n in names}
    ld = pd.read_csv(path, sep="\t").drop_duplicates("name", keep="last").set_index("name")
    return {n: (bool(ld.loc[n, "h2_z"] >= 4) if n in ld.index else None) for n in names}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="decodeme_gwas_1")
    ap.add_argument("--panel-exclude", nargs="*", default=[], help="panel traits to leave out of the null")
    args = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)

    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")

    names = traits_with_results("A")
    panel = [n for n in names if not n.startswith("decodeme_") and n not in args.panel_exclude]
    gate = passes_heritability_gate(panel)
    panel = [n for n in panel if gate.get(n) is not False]
    print("target:", args.target, "| panel traits:", len(panel))

    model_a = load_model(args.target, "A")
    out = ann.join(model_a[["NGENES", "BETA", "SE", "P"]].rename(columns={"P": "p_A"}))
    out["fdr_A"] = multipletests(out["p_A"], method="fdr_bh")[1]
    if (DERIVED / "magma" / args.target / "cluster_B.gsa.out").exists():
        out["p_B"] = load_model(args.target, "B")["P"]
    out["z"] = z_from_p(out["p_A"])
    out["z_std"] = robust_standardise(out["z"])

    # profile of every panel trait across the 461 clusters
    profiles = pd.DataFrame({t: robust_standardise(z_from_p(load_model(t, "A")["P"])) for t in panel},
                            index=out.index)
    out["panel_mean"] = profiles.mean(axis=1)
    out["panel_sd"] = profiles.std(axis=1, ddof=1)
    out["s"] = (out["z_std"] - out["panel_mean"]) / out["panel_sd"]
    out["panel_percentile"] = (profiles.lt(out["z_std"], axis=0)).mean(axis=1)

    out["T2_specific"] = (out["fdr_A"] < 0.05) & (out.get("p_B", 1.0) < 0.05) & (out["s"] >= 2.0)
    out.to_csv(RESULTS / f"clusters_{args.target}.tsv", sep="\t", float_format="%.5g")

    top = out.sort_values("p_A").head(15)
    cols = ["cluster_name", "supercluster", "z", "fdr_A", "s", "panel_percentile", "T2_specific"]
    print(top[cols].to_string())
    print("T2-specific clusters:", int(out["T2_specific"].sum()))
    out.groupby("supercluster").agg(n=("z", "size"), n_fdr=("fdr_A", lambda x: int((x < 0.05).sum())),
                                    n_specific=("T2_specific", "sum"),
                                    max_z=("z", "max")).sort_values("max_z", ascending=False
                                    ).to_csv(RESULTS / f"superclusters_{args.target}.tsv", sep="\t", float_format="%.4g")


if __name__ == "__main__":
    main()
