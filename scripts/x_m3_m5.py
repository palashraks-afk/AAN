"""NEXT_20 M5 (cross-atlas identity of Amex_153/175) and M3 (replicated-loci genes), PREREG v1.15.

    python scripts/x_m3_m5.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DERIVED = Path("D:/AAN_data/derived")
RESULTS = Path("D:/AAN/results")
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]
GENES = {"CLYBL": 171425, "BICD1": 636, "GRIN2A": 2903, "CSMD1": 64478, "RORA": 6095}  # Entrez ids


def main():
    sil = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t").set_index("GENE")
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_name")
    ids = {n: int(ann.loc[n, "cluster_id"]) if "cluster_id" in ann.columns else int(ann.loc[n].name) for n in SIX}
    rows = []
    for name in ("tran", "cao"):
        cov = pd.read_csv(DERIVED / f"atlas_{name}_covar.txt", sep="\t").set_index("GENE")
        names = pd.read_csv(DERIVED / f"atlas_{name}_names.tsv", sep="\t", header=None, index_col=0)[1]
        cols = [c for c in cov.columns if c != "avg_all"]
        shared = sil.index.intersection(cov.index)
        for cl in ("Amex_153", "Amex_175"):
            x = sil.loc[shared, f"c{ids[cl]}"].to_numpy()
            res = sorted(((spearmanr(x, cov.loc[shared, c].to_numpy())[0], names[c]) for c in cols), reverse=True)
            for rank, (r, g) in enumerate(res[:5], 1):
                rows.append({"atlas": name, "siletti_cluster": cl, "rank": rank, "best_match": g, "spearman": round(r, 3), "shared_genes": len(shared)})
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "x_M5_cross_atlas.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    # M3: replicated-loci genes vs expression-decile-matched random sets
    rng = np.random.default_rng(20261010)
    dec = pd.qcut(sil["avg_all"].rank(method="first"), 10, labels=False)
    mine = [g for g in GENES.values() if g in sil.index]
    print("replicated-loci genes present in the specificity matrix:", len(mine), "of", len(GENES))
    pools = {d: sil.index[dec == d].to_numpy() for d in range(10)}
    res = []
    for cl in ("Amex_153", "Amex_175"):
        col = f"c{ids[cl]}"
        obs = sil.loc[mine, col].mean()
        null = np.array([np.mean([sil.loc[rng.choice(pools[dec.loc[g]]), col] for g in mine]) for _ in range(10000)])
        p = (1 + (null >= obs).sum()) / 10001
        res.append({"cluster": cl, "mean_specificity_of_5_genes": round(obs, 3), "null_95th_pct": round(np.percentile(null, 95), 3), "p_one_sided": round(p, 4), "supported_p_lt_0.025": bool(p < 0.025)})
    r = pd.DataFrame(res)
    r.to_csv(RESULTS / "x_M3_replicated_loci_genes.tsv", sep="\t", index=False)
    print(r.to_string(index=False))


if __name__ == "__main__":
    main()
