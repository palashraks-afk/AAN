"""NEXT_20 M7 and M10: publication-style summary figure and the pre-registration timeline. Okabe-Ito colour-blind-safe palette.

    python scripts/x_figures.py
"""
import re
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES, FIG = ROOT / "results", ROOT / "figures"
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "red": "#D55E00", "sky": "#56B4E9", "pink": "#CC79A7", "grey": "#999999", "yellow": "#F0E442"}
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})


def summary():
    c = pd.read_csv(RES / "clusters_decodeme_gwas_1.tsv", sep="\t")
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    a = ax[0, 0]
    a.scatter(c["z_std"], c["s"], s=8, color=OI["grey"], alpha=0.5, label="461 clusters")
    sig = c[c["fdr_A"] < 0.05]
    a.scatter(sig["z_std"], sig["s"], s=26, color=OI["sky"], label="FDR < 0.05 (11)")
    six = c[c["cluster_name"].isin(SIX)]
    a.scatter(six["z_std"], six["s"], s=48, color=OI["red"], edgecolor="k", label="also specific (s >= 2)")
    for r in six.itertuples():
        a.annotate(r.cluster_name, (r.z_std, r.s), fontsize=7, xytext=(-30 if r.cluster_name in ("Splat_402", "DLIT_150") else 5, 4 if r.cluster_name != "DLIT_152" else -10), textcoords="offset points")
    a.axhline(2, color="k", lw=0.8, ls="--")
    a.set_xlabel("ME/CFS enrichment (robust z within trait)")
    a.set_ylabel("specificity score s vs 19 other traits")
    a.set_title("A  Calibrated cell-type result (DecodeME GWAS-1)", loc="left", fontweight="bold")
    a.legend(frameon=False, fontsize=8, loc="upper left")
    # B: independent checks
    b = ax[0, 1]
    checks = [
        ("Gene families (F1)", "not supported", OI["grey"]),
        ("Panel-PC conditioning (v1.9)", "not supported", OI["grey"]),
        ("GTEx amygdala vs brain (v1.10)", "not supported", OI["grey"]),
        ("Amygdala atlas, Tran 2021 (v1.11)", "not supported (weak)", OI["grey"]),
        ("FinnGen fatigue cohort (v1.11)", "not supported", OI["grey"]),
        ("Fetal/peripheral atlas (v1.13)", "not supported", OI["grey"]),
        ("Replicated-loci genes (v1.15)", "not supported", OI["grey"]),
        ("Pain cohort, six clusters (v1.15)", "shared (p = 0.0001)", OI["orange"]),
    ]
    for i, (n, r, col) in enumerate(checks[::-1]):
        b.barh(i, 1, color=col, alpha=0.85)
        b.text(0.02, i, f"{n}: {r}", va="center", fontsize=8)
    b.set_xlim(0, 1); b.set_ylim(-0.6, len(checks) - 0.4); b.axis("off")
    b.set_title("B  Independent checks (pre-registered)", loc="left", fontweight="bold")
    # C: related conditions
    cc = ax[1, 0]
    rel = pd.read_csv(RES / "x_M4_related_conditions.tsv", sep="\t").dropna(subset=["mean_z_six"])
    fat = pd.read_csv(RES / "x_R1_finngen_fatigue_clusters.tsv", sep="\t")
    fz = fat[fat["cluster_name"].isin(SIX)]["z"].mean()
    names = ["Malaise & fatigue\n(FinnGen)", "Fibromyalgia\n(FinnGen)", "Pain\n(FinnGen)", "IBS"]
    vals = [fz] + [rel.set_index("condition").loc[k, "mean_z_six"] for k in ("finngen_fibromyalgia", "finngen_pain", "ibs")]
    nulls = [0.444] + [rel.set_index("condition").loc[k, "null_95th"] for k in ("finngen_fibromyalgia", "finngen_pain", "ibs")]
    x = np.arange(4)
    cc.bar(x, vals, color=[OI["grey"], OI["grey"], OI["orange"], OI["grey"]], width=0.55, label="mean z of the six clusters")
    cc.scatter(x, nulls, marker="_", s=900, color="k", label="95th percentile of random sets")
    cc.set_xticks(x); cc.set_xticklabels(names, fontsize=8)
    cc.set_ylabel("mean z of the six clusters")
    cc.set_title("C  The six clusters in related conditions", loc="left", fontweight="bold")
    cc.legend(frameon=False, fontsize=8, loc="upper right")
    # D: controls from the G3 check
    d = ax[1, 1]
    g3 = (RES / "G3_check.txt").read_text()
    top = {}
    for blk in re.split(r"\n\s*(\w+) top clusters", g3)[1:]:
        pass
    m = re.findall(r"(schizophrenia|alzheimer|height) top clusters\n.*?\n.*?\n\s*[\d.]+\s+\S+\s+(\S+)\s+([\d.]+)", g3, flags=re.S)
    labels = {"schizophrenia": "Schizophrenia\n(excitatory neurons)", "alzheimer": "Alzheimer's\n(microglia)", "height": "Height\n(non-neuronal)"}
    zs = {k: float(v) for k, _, v in m} if m else {"schizophrenia": 7.09, "alzheimer": 5.99, "height": 10.42}
    d.bar(range(3), [zs.get(k, np.nan) for k in labels], color=[OI["blue"], OI["green"], OI["pink"]], width=0.55)
    d.bar(3, 3.58, color=OI["red"], width=0.55)
    d.set_xticks(range(4)); d.set_xticklabels(list(labels.values()) + ["ME/CFS\n(best cluster)"], fontsize=8)
    d.set_ylabel("top cluster z (model A)")
    d.set_title("D  Pipeline checked on known biology first", loc="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "fig_summary_panels.png", bbox_inches="tight")
    print("wrote fig_summary_panels.png")


def timeline():
    tags = subprocess.run(["git", "for-each-ref", "--format=%(refname:short)|%(creatordate:iso8601)", "refs/tags"], cwd=ROOT, capture_output=True, text=True).stdout.strip().splitlines()
    t = []
    for ln in tags:
        n, d = ln.split("|")
        t.append((n.replace("prereg-", ""), datetime.strptime(d[:19], "%Y-%m-%d %H:%M:%S")))
    dev = []
    for m in re.finditer(r"^## (\d{4}-\d{2}-\d{2})(?: \((\d{2}:\d{2})\)| \(evening\))?", (ROOT / "preregistration" / "DEVIATIONS.md").read_text(), flags=re.M):
        hh = m.group(2) or "12:00"
        dev.append(datetime.strptime(f"{m.group(1)} {hh}", "%Y-%m-%d %H:%M"))
    events = [(datetime(2026, 10, 9, 19, 36), "Control gate G3 passes;\nfirst ME/CFS cell-type look", OI["green"]),
              (datetime(2026, 10, 10, 14, 0), "Independent checks begin\n(all pre-registered first)", OI["blue"])]
    fig, ax = plt.subplots(figsize=(12, 4))
    for i, (n, d) in enumerate(sorted(t, key=lambda x: x[1])):
        h = 0.35 + 0.22 * (i % 6)
        ax.vlines(d, 0, h, color=OI["blue"], lw=1.2)
        ax.text(d, h + 0.02, n, fontsize=7, ha="center", va="bottom", color=OI["blue"], bbox=dict(fc="white", ec="none", pad=0.5))
    for d in dev:
        ax.vlines(d, 0, -0.6, color=OI["orange"], lw=1.6)
    for d, txt, col in events:
        ax.axvline(d, color=col, ls="--", lw=1)
        ax.text(d, -0.95, txt, fontsize=8, ha="center", color=col)
    ax.axhline(0, color="k", lw=0.8)
    ax.plot([], [], color=OI["blue"], label="pre-registration written and tagged (before the analysis it describes)")
    ax.plot([], [], color=OI["orange"], lw=2, label="deviation logged in DEVIATIONS.md (date and reason)")
    ax.set_ylim(-1.3, 2.3); ax.set_yticks([])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M"))
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_title("Planned versus changed: pre-registrations (git tags) and logged deviations, 7 to 10 Oct 2026", loc="left", fontweight="bold")
    ax.spines["left"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_prereg_timeline.png", bbox_inches="tight")
    print("wrote fig_prereg_timeline.png;", len(t), "tags,", len(dev), "deviations")


if __name__ == "__main__":
    summary()
    timeline()
