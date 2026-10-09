"""Discovery layer (D1 to D4) and convergence/reporting layers (E1, E2), as defined in preregistration v1.2 and v1.3.

Run only after the control check (check_controls.py) has passed.

    python scripts/analyze_layers.py
Writes tables to results/ and quantile-quantile plots to figures/.
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr
from statsmodels.stats.multitest import multipletests

import brain3d
import calibrate

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DERIVED = Path("D:/AAN_data/derived")
RES = Path("D:/AAN/results")
FIG = Path("D:/AAN/figures")
TARGET, ROBUST = "decodeme_gwas_1", "decodeme_gwas_2"
LITERATURE_REGIONS = ["Basal nuclei", "Thalamus", "Midbrain", "Pons", "Medulla", "Hippocampus", "Amygdala",
                      "Hypothalamus"]                          # fixed in preregistration v1.3 (E1)


def gsa(name, prefix, id_pattern=r".+"):
    for suffix in (".gsa.out.txt", ".gsa.out"):
        path = DERIVED / "magma" / name / f"{prefix}{suffix}"
        if path.exists():
            df = pd.read_csv(path, comment="#", sep=r"\s+")
            df = df[df["VARIABLE"].str.match(id_pattern) & ~df["VARIABLE"].isin(["avg_all", "avg_neuron"])]
            return df.set_index("VARIABLE")
    return None


def panel_traits(prefix):
    out = []
    for d in sorted((DERIVED / "magma").iterdir()):
        if d.name.startswith("decodeme_"):
            continue
        if gsa(d.name, prefix) is not None:
            out.append(d.name)
    return out


def calibrated(prefix, label):
    """Same calibration as the main analysis, for group or cell-type tables."""
    tgt = gsa(TARGET, prefix)
    rob = gsa(ROBUST, prefix)
    out = tgt[["NGENES", "BETA", "SE", "P"]].copy()
    out["z"] = calibrate.z_from_p(out["P"])
    out["fdr"] = multipletests(out["P"], method="fdr_bh")[1]
    out["z_std"] = calibrate.robust_standardise(out["z"])
    panel = {t: calibrate.robust_standardise(calibrate.z_from_p(gsa(t, prefix)["P"])) for t in panel_traits(prefix)}
    prof = pd.DataFrame(panel, index=out.index)
    out["panel_mean"], out["panel_sd"] = prof.mean(axis=1), prof.std(axis=1, ddof=1)
    out["s"] = (out["z_std"] - out["panel_mean"]) / out["panel_sd"]
    if rob is not None:
        out["p_gwas2"] = rob["P"]
        out["z_gwas2"] = calibrate.z_from_p(rob["P"])
    out["supported"] = (out["fdr"] < 0.05) & (out["s"] >= 2.0) & (out.get("p_gwas2", 1.0) < 0.05)
    out.to_csv(RES / f"{label}_{TARGET}.tsv", sep="\t", float_format="%.5g")
    print(f"\n{label}: {int((out['fdr'] < 0.05).sum())} of {len(out)} FDR<0.05, {int(out['supported'].sum())} supported")
    print(out.sort_values("P").head(8)[["z", "fdr", "s", "supported"]].to_string())
    return out


def pathways():
    names = pd.read_csv(DERIVED / "reactome_names.tsv", sep="\t").set_index("set")
    tables = {}
    for prefix in ("reactome_A", "reactome_C"):
        for who in (TARGET, ROBUST):
            t = gsa(who, prefix)
            if t is not None:
                t["fdr"] = multipletests(t["P"], method="fdr_bh")[1]
                tables[(who, prefix)] = t
    a, c = tables.get((TARGET, "reactome_A")), tables.get((TARGET, "reactome_C"))
    if a is None:
        return
    out = a[["NGENES", "BETA", "P", "fdr"]].rename(columns={"P": "p_A", "fdr": "fdr_A"}).join(names["title"])
    if c is not None:
        out = out.join(c[["P", "fdr"]].rename(columns={"P": "p_C", "fdr": "fdr_C"}))
    r = tables.get((ROBUST, "reactome_A"))
    if r is not None:
        out["p_gwas2"] = r["P"]
    out.sort_values("p_A").to_csv(RES / f"pathways_{TARGET}.tsv", sep="\t", float_format="%.5g")
    print("\npathways FDR<0.05 (model A):", int((out["fdr_A"] < 0.05).sum()))
    print(out.sort_values("p_A").head(10)[["title", "p_A", "fdr_A"]].to_string())


def comorbidity_map(clusters):
    traits = calibrate.traits_with_results("A")
    z = {t: calibrate.z_from_p(calibrate.load_model(t, "A")["P"]) for t in traits if not t.startswith("decodeme_gwas_1_")
         and t != "decodeme_gwas_2"}
    z["ME/CFS"] = clusters["z"]
    mat = pd.DataFrame(z)
    corr = mat.corr(method="spearman")
    corr.to_csv(RES / "profile_correlation_matrix.tsv", sep="\t", float_format="%.4f")
    others = [c for c in corr.columns if c != "ME/CFS"]
    pair_vals = [corr.loc[a, b] for i, a in enumerate(others) for b in others[i + 1:]]
    me = corr.loc["ME/CFS", others].sort_values(ascending=False)
    ref = pd.Series(pair_vals)
    me_df = pd.DataFrame({"spearman_with_mecfs": me,
                          "percentile_among_other_trait_pairs": [(ref < v).mean() for v in me]})
    me_df.to_csv(RES / "comorbidity_map.tsv", sep="\t", float_format="%.4f")
    print("\ncellular comorbidity map (top 8 by correlation with ME/CFS):")
    print(me_df.head(8).to_string())


def convergence(clusters):
    def stat(scores):
        inside = scores.reindex(LITERATURE_REGIONS).dropna()
        outside = scores.drop(inside.index)
        return inside.mean() - outside.mean()

    me_scores = brain3d.region_scores(clusters["z"])
    me = stat(me_scores)
    traits = [t for t in calibrate.traits_with_results("A") if not t.startswith("decodeme_")]
    panel = {t: stat(brain3d.region_scores(calibrate.z_from_p(calibrate.load_model(t, "A")["P"]))) for t in traits}
    pct = float(np.mean([v < me for v in panel.values()]))
    rng = np.random.default_rng(20261008)
    k = len(me_scores.reindex(LITERATURE_REGIONS).dropna())
    perm = [me_scores.iloc[rng.permutation(len(me_scores))[:k]].mean()
            - me_scores.iloc[rng.permutation(len(me_scores))[k:]].mean() for _ in range(10_000)]
    p_perm = float((np.sum(np.array(perm) >= me) + 1) / (len(perm) + 1))
    pd.DataFrame({"trait": ["ME/CFS"] + list(panel), "inside_minus_outside": [me] + list(panel.values())}).to_csv(
        RES / "convergence_E1.tsv", sep="\t", index=False, float_format="%.4f")
    verdict = "converges" if (pct >= 0.90 and p_perm < 0.05) else "no convergence beyond generic brain traits"
    print(f"\nE1 literature-region convergence: statistic {me:.3f}; percentile among panel {pct:.2f}; "
          f"permutation p {p_perm:.4f} -> {verdict}")


def qq_and_lambda(clusters):
    genes = pd.read_csv(DERIVED / "magma" / TARGET / "genes.genes.out", sep="\t")
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6))
    for ax, (label, p) in zip(axes, (("gene-level p (MAGMA)", genes["P"].to_numpy()),
                                      ("cluster-level p (model A)", clusters["P"].to_numpy()))):
        p = np.sort(p[p > 0])
        obs, exp = -np.log10(p), -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
        ax.scatter(exp, obs, s=4, c="#34495e", linewidth=0)
        ax.plot([0, exp.max()], [0, exp.max()], "k", lw=0.6)
        ax.set_xlabel("expected -log10 p"); ax.set_ylabel("observed -log10 p"); ax.set_title(label, fontsize=9)
        lam = float(np.median(norm.isf(p / 2) ** 2) / 0.4549)
        rows.append({"what": label, "lambda_GC": lam, "n": len(p)})
    fig.savefig(FIG / "figS_qq.png", dpi=300, bbox_inches="tight")
    pd.DataFrame(rows).to_csv(RES / "inflation_E2.tsv", sep="\t", index=False, float_format="%.3f")
    print("\nE2 inflation:", rows)


def main():
    RES.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    clusters = pd.read_csv(RES / f"clusters_{TARGET}.tsv", sep="\t", index_col="cluster_id")
    calibrated("groups_A", "groups_A")
    calibrated("groups_B", "groups_B")
    calibrated("hpa_A", "hpa")
    pathways()
    comorbidity_map(clusters)
    convergence(clusters)
    qq_and_lambda(clusters)


if __name__ == "__main__":
    main()
