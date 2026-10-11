"""Figures for the final report (all read from files in results/ and D:/AAN_data/derived)."""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

import calibrate  # noqa: E402

RES = Path("D:/AAN/results")
FIG = Path("D:/AAN/figures")
DER = Path("D:/AAN_data/derived")
INK, GREY, ACC, ACC2, BLUE = "#2c3e50", "#95a5a6", "#c0392b", "#e67e22", "#2874a6"
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]
ORDER = ["Amex_153", "Amex_175", "ULIT_121", "Splat_402", "DLIT_152", "DLIT_150"]
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})


def pipeline():
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.axis("off")
    boxes = [("DecodeME GWAS\n15,579 cases\n259,909 controls", 0.0), ("SNP p-values\nto gene scores\n(MAGMA)", 0.2),
             ("Cell-type test\n461 brain clusters\n(Siletti atlas)", 0.4), ("Compare with\n19 other traits\n(specificity score s)", 0.6),
             ("Checks: permutation,\nlocus drop, 2nd GWAS,\nconditioning, windows", 0.8)]
    for text, x in boxes:
        ax.add_patch(FancyBboxPatch((x, 0.4), 0.17, 0.4, boxstyle="round,pad=0.01", fc="#eaf2f8", ec=BLUE, lw=1.5, transform=ax.transAxes))
        ax.text(x + 0.085, 0.6, text, ha="center", va="center", fontsize=10, color=INK, transform=ax.transAxes)
    for x in (0.18, 0.38, 0.58, 0.78):
        ax.annotate("", xy=(x + 0.02, 0.6), xytext=(x, 0.6), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color=INK, lw=1.5))
    ax.add_patch(FancyBboxPatch((0.02, 0.02), 0.93, 0.22, boxstyle="round,pad=0.01", fc="#fdf2e9", ec=ACC2, lw=1.5, transform=ax.transAxes))
    ax.text(0.485, 0.13, "Before any ME/CFS cell-type result was read, control traits had to reproduce known biology\n"
            "(schizophrenia: neurons, Alzheimer's: microglia, height: connective tissue). The plan was written and tagged in advance.",
            ha="center", va="center", fontsize=9, color=INK, transform=ax.transAxes)
    fig.savefig(FIG / "fig0_pipeline.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def calibration():
    c = pd.read_csv(RES / "clusters_decodeme_gwas_1.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(c["panel_mean"], c["z_std"], s=9, c=GREY, linewidth=0, alpha=0.8, label="other clusters")
    sig = c[c["fdr_A"] < 0.05]
    ax.scatter(sig["panel_mean"], sig["z_std"], s=28, c=BLUE, label="enriched in ME/CFS (FDR < 0.05)")
    sp = c[c["T2_specific"]]
    ax.scatter(sp["panel_mean"], sp["z_std"], s=60, c=ACC, edgecolor="white", linewidth=0.8, label="also specific (s >= 2)")
    offs = {"Amex_153": (-95, 30), "Amex_175": (-120, 12), "ULIT_121": (-30, -42), "Splat_402": (10, -45)}
    for _, r in sp.iterrows():
        if r["cluster_name"] in offs:
            ax.annotate(r["cluster_name"], (r["panel_mean"], r["z_std"]), xytext=offs[r["cluster_name"]], textcoords="offset points",
                        fontsize=9, color=INK, arrowprops=dict(arrowstyle="-", color=GREY, lw=0.8))
    lim = [-2.5, 2.6]
    ax.plot(lim, lim, c=INK, lw=0.7, ls="--")
    ax.set_xlabel("average enrichment of the 19 comparison traits (standardised)")
    ax.set_ylabel("ME/CFS enrichment (standardised)")
    ax.set_title("Is the ME/CFS profile different from other traits?", fontsize=12)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.savefig(FIG / "fig3_calibration_clean.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def controls():
    ann = pd.read_csv(DER / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    traits = [("schizophrenia", "Schizophrenia"), ("alzheimer", "Alzheimer's"), ("height", "Height"), ("decodeme_gwas_1", "ME/CFS")]
    rows = {}
    for t, lab in traits:
        m = calibrate.load_model(t, "A")
        z = pd.Series(calibrate.z_from_p(m["P"]), index=m.index)
        rows[lab] = z.groupby(ann.loc[z.index, "supercluster"]).max()
    df = pd.DataFrame(rows).fillna(0)
    keep = df.max(axis=1).sort_values(ascending=False).head(14).index
    df = df.loc[keep]
    fig, ax = plt.subplots(figsize=(7.5, 6))
    im = ax.imshow(df.values, cmap="YlOrRd", aspect="auto", vmin=0)
    ax.set_xticks(range(df.shape[1]), df.columns)
    ax.set_yticks(range(df.shape[0]), df.index)
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            ax.text(j, i, f"{df.values[i, j]:.1f}", ha="center", va="center", fontsize=8,
                    color="black" if df.values[i, j] < 5 else "white")
    ax.set_title("Controls: strongest enrichment (z) per cell class", fontsize=12)
    fig.colorbar(im, ax=ax, label="z", shrink=0.8)
    ax.spines[:].set_visible(False)
    fig.savefig(FIG / "fig2_controls.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def validation():
    d = pd.read_csv(RES / "panel_as_target.tsv", sep="\t")
    d = d[d["n_sig"] >= 3].sort_values("share_specific")
    labels = list(d["trait"]) + ["ME/CFS"]
    vals = list(d["share_specific"] * 100) + [6 / 11 * 100]
    ns = list(d["n_sig"]) + [11]
    cols = [BLUE] * len(d) + [ACC]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.barh(labels, vals, color=cols)
    for i, (v, n) in enumerate(zip(vals, ns)):
        ax.text(v + 1, i, f"{v:.0f}%  ({n} enriched clusters)", va="center", fontsize=9, color=INK)
    ax.set_xlim(0, 115)
    ax.set_xlabel("share of enriched clusters that are specific to the trait (%)")
    ax.set_title("Does the specificity score work? Each trait scored against the others", fontsize=11)
    fig.savefig(FIG / "fig5_validation.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def robustness():
    r = pd.read_csv(RES / "panel_robustness.tsv", sep="\t").set_index("cluster").loc[ORDER]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    series = [("full panel", "full panel", INK, "o"), ("R4 s without similar", "without similar traits", BLUE, "s"),
              ("R4 s without psych", "without psychiatric traits", ACC2, "^"), ("R2 min s", "worst leave-one-trait-out", ACC, "v")]
    for col, lab, c, mk in series:
        ax.scatter(r[col], range(len(r)), c=c, marker=mk, s=45, label=lab, zorder=3)
    ax.axvline(2, c=GREY, ls="--", lw=1)
    ax.set_yticks(range(len(r)), r.index)
    ax.invert_yaxis()
    ax.set_xlabel("specificity score s (2 = the cutoff)")
    ax.set_title("Specificity under different comparison panels", fontsize=12)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.savefig(FIG / "fig6_robustness.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def subsets():
    ann = pd.read_csv(DER / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    ids = ann.index[ann["cluster_name"].isin(SIX)]
    cols = {"all cases": None, "infectious onset": "infectious_onset", "non-infectious onset": "non_infectious_onset",
            "female": "female", "male": "male"}
    M = {}
    for lab, s in cols.items():
        t = "decodeme_gwas_1" if s is None else f"decodeme_gwas_1_{s}"
        m = calibrate.load_model(t, "A")
        M[lab] = -np.log10(m.loc[ids, "P"].values)
    df = pd.DataFrame(M, index=ann.loc[ids, "cluster_name"].values).loc[ORDER]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    im = ax.imshow(df.values, cmap="YlOrRd", aspect="auto", vmin=0, vmax=4.5)
    ax.set_xticks(range(df.shape[1]), df.columns, rotation=20, ha="right")
    ax.set_yticks(range(df.shape[0]), df.index)
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            ax.text(j, i, f"{df.values[i, j]:.1f}", ha="center", va="center", fontsize=9,
                    color="black" if df.values[i, j] < 3 else "white")
    ax.set_title("Signal in the six clusters by patient subset (-log10 p)", fontsize=11)
    fig.colorbar(im, ax=ax, shrink=0.8, label="-log10 p")
    ax.spines[:].set_visible(False)
    fig.savefig(FIG / "fig9_subsets.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def windows_and_genes():
    ann = pd.read_csv(DER / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    ids = ann.index[ann["cluster_name"].isin(SIX)]
    res = {}
    for lab, name in (("35/10 kb (planned)", "decodeme_gwas_1"), ("0 kb", "decodeme_gwas_1_w0"), ("50 kb", "decodeme_gwas_1_w50")):
        m = calibrate.load_model(name, "A")
        res[lab] = pd.Series(calibrate.z_from_p(m["P"]), index=m.index).loc[ids].values
    w = pd.DataFrame(res, index=ann.loc[ids, "cluster_name"].values).loc[ORDER]
    g = pd.read_csv(RES / "w9_drivers.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.arange(len(w))
    for k, (col, c) in enumerate(zip(w.columns, (INK, BLUE, ACC2))):
        axes[0].bar(x + (k - 1) * 0.27, w[col], 0.27, label=col, color=c)
    axes[0].set_xticks(x, w.index, rotation=25, ha="right")
    axes[0].set_ylabel("enrichment z")
    axes[0].set_title("Choice of gene window", fontsize=11)
    axes[0].set_ylim(0, 4.4)
    axes[0].legend(frameon=False, fontsize=8.5, ncol=3, loc="upper center")
    for name, c in (("Amex_153", ACC), ("Amex_175", BLUE)):
        s = g[g.cluster == name]
        axes[1].plot(s["genes_dropped"], s["z"], marker="o", c=c, label=name)
    zp = float(calibrate.z_from_p(0.05))
    axes[1].axhline(zp, c=GREY, ls="--", lw=1)
    axes[1].text(0.3, zp + 0.05, "p = 0.05", color=GREY, fontsize=8.5)
    axes[1].set_xlabel("top contributing genes removed")
    axes[1].set_ylabel("enrichment z")
    axes[1].set_title("Does the signal rest on a few genes?", fontsize=11)
    axes[1].legend(frameon=False)
    fig.savefig(FIG / "fig10_windows_genes.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    for f in (pipeline, calibration, controls, validation, robustness, subsets, windows_and_genes):
        f()
        print("done", f.__name__)
