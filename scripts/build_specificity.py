"""Turn the aggregated Human Brain Cell Atlas into gene-by-cluster specificity scores for MAGMA.

Input : adult_human_agg.loom (59,480 genes x 461 clusters, mean UMI per cell)
        cluster_annotation.xlsx, NCBI38.gene.loc
Output: derived/gene_covar.txt      Entrez id + one column per cluster + two control columns
        derived/cluster_annotation.tsv
        derived/top10_sets.txt      top 10% specific genes per cluster (sensitivity check)
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

# superclusters that are not neurons; everything else counts as neuronal
NON_NEURONAL = {
    "Astrocyte", "Microglia", "Oligodendrocyte", "Oligodendrocyte precursor",
    "Committed oligodendrocyte precursor", "Vascular", "Fibroblast", "Ependymal",
    "Choroid plexus", "Bergmann glia", "Miscellaneous",
}
MIN_MEAN_CPM = 1.0   # genes below this average across clusters are dropped, specificity is noise there


def rank_normal(x):
    r = rankdata(x, method="average")
    return norm.ppf((r - 0.5) / len(x))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default="D:/AAN_data/atlas/adult_human_agg.loom")
    ap.add_argument("--annot", default="D:/AAN_data/atlas/cluster_annotation.xlsx")
    ap.add_argument("--geneloc", default="D:/AAN_data/tools/ref/NCBI38.gene.loc")
    ap.add_argument("--outdir", default="D:/AAN_data/derived")
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.atlas, "r") as f:
        mat = f["matrix"][:]
        symbols = np.array([g.decode() for g in f["row_attrs"]["Gene"][:]])
        cluster_ids = f["col_attrs"]["Clusters"][:]

    ann = pd.read_excel(args.annot)
    ann.columns = [c.strip() for c in ann.columns]
    ann = ann.rename(columns={"Cluster ID": "cluster_id", "Cluster name": "cluster_name",
                              "Supercluster": "supercluster"})
    ann = ann.set_index("cluster_id").loc[cluster_ids].reset_index()
    ann["neuronal"] = ~ann["supercluster"].isin(NON_NEURONAL)
    ann[["cluster_id", "cluster_name", "supercluster", "neuronal", "Number of cells"]].to_csv(
        out / "cluster_annotation.tsv", sep="\t", index=False)

    # counts per million within each cluster, so depth differences between clusters don't matter
    cpm = mat / mat.sum(axis=0, keepdims=True) * 1e6

    loc = pd.read_csv(args.geneloc, sep="\t", header=None,
                      names=["entrez", "chr", "start", "end", "strand", "symbol"])
    loc = loc.drop_duplicates("symbol")
    sym_to_entrez = dict(zip(loc["symbol"], loc["entrez"]))

    # a symbol can appear twice in the atlas, keep the row with the most expression
    df = pd.DataFrame(cpm, columns=cluster_ids)
    df["symbol"] = symbols
    df["total"] = cpm.sum(axis=1)
    df = df.sort_values("total", ascending=False).drop_duplicates("symbol")
    df = df[df["symbol"].isin(sym_to_entrez)]
    df = df[df["total"] / len(cluster_ids) >= MIN_MEAN_CPM]
    print("genes kept:", len(df))

    expr = df[list(cluster_ids)].to_numpy()
    entrez = df["symbol"].map(sym_to_entrez).to_numpy()

    spec = expr / expr.sum(axis=1, keepdims=True)           # share of the gene's expression in each cluster
    spec_rn = np.apply_along_axis(rank_normal, 0, spec)     # rank-normal within each cluster

    neuronal = ann["neuronal"].to_numpy()
    avg_all = np.log10(expr.mean(axis=1) + 1)
    avg_neuron = np.log10(expr[:, neuronal].mean(axis=1) + 1)

    covar = pd.DataFrame(spec_rn, columns=[f"c{c}" for c in cluster_ids])
    covar.insert(0, "GENE", entrez)
    covar["avg_all"] = avg_all
    covar["avg_neuron"] = avg_neuron
    covar.to_csv(out / "gene_covar.txt", sep="\t", index=False, float_format="%.5f")

    # top 10% of genes by specificity in each cluster, for the sensitivity analysis
    cutoff = int(round(0.10 * len(df)))
    with open(out / "top10_sets.txt", "w") as fh:
        for j, cid in enumerate(cluster_ids):
            top = np.argsort(spec[:, j])[::-1][:cutoff]
            fh.write(f"c{cid}\t" + " ".join(str(g) for g in entrez[top]) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
