"""L2 input (PREREG v1.5): gene specificity across the 17 brain regions of the atlas.

Region expression is the cell-count-weighted mean of the cluster expression of the clusters found in that region's
dissections. Output: gene_covar_regions.txt (same genes and the avg_all column as the cluster file).
"""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

import brain3d

ATLAS = Path("D:/AAN_data/atlas/adult_human_agg.loom")
DERIVED = Path("D:/AAN_data/derived")
GENELOC = "D:/AAN_data/tools/ref/NCBI38.gene.loc"


def rank_normal(x):
    r = rankdata(x, method="average")
    return norm.ppf((r - 0.5) / len(x))


def main():
    covar = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t", usecols=["GENE", "avg_all"])
    counts = brain3d.cluster_by_region_counts()                       # clusters x regions
    with h5py.File(ATLAS, "r") as f:
        mat = f["matrix"][:]
        symbols = np.array([g.decode() for g in f["row_attrs"]["Gene"][:]])
        ids = list(f["col_attrs"]["Clusters"][:])
    cpm = pd.DataFrame(mat / mat.sum(axis=0, keepdims=True) * 1e6, columns=ids)
    loc = pd.read_csv(GENELOC, sep="\t", header=None, names=["entrez", "chr", "s", "e", "st", "symbol"])
    sym2ent = dict(zip(loc["symbol"].drop_duplicates(), loc.drop_duplicates("symbol")["entrez"]))
    cpm["entrez"] = [sym2ent.get(s) for s in symbols]
    cpm = cpm.dropna(subset=["entrez"])
    cpm = cpm.assign(total=cpm[ids].sum(axis=1)).sort_values("total", ascending=False).drop_duplicates("entrez")
    cpm = cpm.set_index(cpm["entrez"].astype(int))[ids].loc[covar["GENE"].values]

    w = counts.loc[ids].to_numpy(dtype=float)                           # clusters x regions
    region_expr = cpm.to_numpy() @ w / w.sum(axis=0, keepdims=True)    # genes x regions
    total = region_expr.sum(axis=1, keepdims=True)
    spec = np.divide(region_expr, total, out=np.zeros_like(region_expr), where=total > 0)
    out = pd.DataFrame(np.apply_along_axis(rank_normal, 0, spec),
                       columns=[c.replace(" ", "_") for c in counts.columns])
    out.insert(0, "GENE", covar["GENE"].values)
    out["avg_all"] = covar["avg_all"].values
    out.to_csv(DERIVED / "gene_covar_regions.txt", sep="\t", index=False, float_format="%.5f")
    print("regions:", list(counts.columns), "genes:", len(out))


if __name__ == "__main__":
    main()
