"""Figures for the report. Run after calibrate.py has written results/clusters_*.tsv.

    python scripts/make_figures.py
Figures land in figures/ as 300 dpi PNG; missing inputs are skipped with a message.
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import calibrate  # noqa: E402

ROOT = Path("D:/AAN")
RES = ROOT / "results"
FIG = ROOT / "figures"
TARGET = "decodeme_gwas_1"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def supercluster_colors(names):
    cmap = plt.get_cmap("tab20")
    return {n: cmap(i % 20) for i, n in enumerate(sorted(set(names)))}


def fig_enrichment(df):
    d = df.sort_values(["supercluster", "cluster_name"]).reset_index()
    cols = supercluster_colors(d["supercluster"])
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.scatter(np.arange(len(d)), d["z"], c=[cols[s] for s in d["supercluster"]], s=10, linewidth=0)
    z_fdr = d.loc[d["fdr_A"] < 0.05, "z"].min() if (d["fdr_A"] < 0.05).any() else np.nan
    if not np.isnan(z_fdr):
        ax.axhline(z_fdr, color="k", lw=0.7, ls="--")
    for _, r in d.nsmallest(6, "p_A").iterrows():
        ax.annotate(r["cluster_name"], (r.name, r["z"]), fontsize=6, xytext=(2, 3), textcoords="offset points")
    centres = d.groupby("supercluster").apply(lambda g: g.index.mean())
    ax.set_xticks(centres.values)
    ax.set_xticklabels(centres.index, rotation=75, ha="right", fontsize=6)
    ax.set_ylabel("enrichment z (MAGMA, model A)")
    ax.set_title("ME/CFS genetic signal across 461 human brain cell clusters")
    fig.savefig(FIG / "fig2_enrichment.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_calibration(df):
    fig, ax = plt.subplots(figsize=(5.6, 5))
    spec = df["T2_specific"]
    ax.fill_between([df["panel_mean"].min() - 1, df["panel_mean"].max() + 1],
                    [df["panel_mean"].min() - 1 - 2 * df["panel_sd"].median()] * 2,
                    [df["panel_mean"].max() + 1 + 2 * df["panel_sd"].median()] * 2, color="#eeeeee", zorder=0)
    ax.scatter(df["panel_mean"], df["z_std"], s=9, c=np.where(spec, "#c0392b", "#7f8c8d"), linewidth=0)
    lo, hi = df[["panel_mean", "z_std"]].min().min(), df[["panel_mean", "z_std"]].max().max()
    ax.plot([lo, hi], [lo, hi], color="k", lw=0.6)
    for _, r in df.nlargest(5, "s").iterrows():
        ax.annotate(r["cluster_name"], (r["panel_mean"], r["z_std"]), fontsize=6, xytext=(3, 2),
                    textcoords="offset points")
    ax.set_xlabel("average standardised enrichment of the comparison panel")
    ax.set_ylabel("ME/CFS standardised enrichment")
    ax.set_title("Is the ME/CFS profile different from other brain traits?")
    fig.savefig(FIG / "fig3_calibration.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_profile_similarity(df):
    names = calibrate.traits_with_results("A")
    panel = calibrate.panel_names(names)
    z = {n: calibrate.z_from_p(calibrate.load_model(n, "A")["P"]) for n in panel}
    z["ME/CFS"] = df["z"]
    mat = pd.DataFrame(z)
    corr = mat.corr(method="spearman")["ME/CFS"].drop("ME/CFS").sort_values()
    fig, ax = plt.subplots(figsize=(4.6, 5.2))
    ax.barh(corr.index, corr.values, color="#34495e")
    ax.set_xlabel("Spearman correlation of the 461-cluster profile with ME/CFS")
    ax.set_title("Which traits share ME/CFS's cell-type profile?")
    fig.savefig(FIG / "fig3b_profile_similarity.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    corr.to_csv(RES / "profile_similarity.tsv", sep="\t", header=["spearman_with_mecfs"])


def fig_conditioning(df):
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    ax.scatter(-np.log10(df["p_A"]), -np.log10(df["p_B"]), s=8, c="#2c3e50", linewidth=0)
    m = max(-np.log10(df["p_A"]).max(), -np.log10(df["p_B"]).max())
    ax.plot([0, m], [0, m], color="k", lw=0.6)
    ax.set_xlabel("-log10 p, model A (all-cluster expression)")
    ax.set_ylabel("-log10 p, model B (+ neuronal expression)")
    ax.set_title("Effect of conditioning on generic neuronal expression")
    fig.savefig(FIG / "fig4_conditioning.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_brain3d(df):
    import brain3d
    scores = brain3d.region_scores(df["z"])
    scores.to_csv(RES / "region_scores_mecfs.tsv", sep="\t", header=["mean_z"])
    brain3d.static_views(scores.drop("Cerebral cortex", errors="ignore"), FIG / "fig7_brain3d.png",
                         "ME/CFS genetic enrichment by brain region (schematic, approximate region centres)")
    brain3d.interactive(scores.drop("Cerebral cortex", errors="ignore"), FIG / "fig7_brain3d.html",
                        "ME/CFS genetic enrichment by brain region")


def main():
    FIG.mkdir(exist_ok=True)
    path = RES / f"clusters_{TARGET}.tsv"
    if not path.exists():
        print("run calibrate.py first:", path)
        return
    df = pd.read_csv(path, sep="\t", index_col="cluster_id")
    for fn in (fig_enrichment, fig_calibration, fig_profile_similarity, fig_conditioning, fig_brain3d):
        try:
            fn(df)
            print("made", fn.__name__)
        except Exception as err:      # keep going so one broken figure doesn't hide the others
            print("FAILED", fn.__name__, type(err).__name__, err)


if __name__ == "__main__":
    main()
