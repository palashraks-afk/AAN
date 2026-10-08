"""Rank genes for one trait and test whether known drug targets rank high (X1 and X3 in the pre-registration).

Gene score = max over the trait's enriched clusters of (MAGMA gene z) x (cluster specificity of the gene).
For a benchmark disease the score is compared with matched random gene sets.

    python scripts/rank_targets.py --trait migraine --disease migraine
    python scripts/rank_targets.py --trait decodeme_gwas_1          # writes the ME/CFS candidate list
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import calibrate

DERIVED = Path("D:/AAN_data/derived")
RESULTS = Path("D:/AAN/results")
N_PERM = 10_000
SEED = 20261007


def gene_table(trait):
    genes = pd.read_csv(DERIVED / "magma" / trait / "genes.genes.out", sep=r"\s+")
    loc = pd.read_csv("D:/AAN_data/tools/ref/NCBI38.gene.loc", sep="\t", header=None,
                      names=["GENE", "chr", "start", "end", "strand", "symbol"])
    return genes.merge(loc[["GENE", "symbol"]], on="GENE", how="left")


def enriched_clusters(trait, min_n=3):
    m = calibrate.load_model(trait, "A")
    fdr = multipletests(m["P"], method="fdr_bh")[1]
    chosen = m.index[fdr < 0.05].tolist()
    if len(chosen) < min_n:      # nothing (or almost nothing) passes: fall back to the strongest few
        chosen = m["P"].nsmallest(min_n).index.tolist()
    return chosen


def score_genes(trait):
    covar = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t")
    genes = gene_table(trait).merge(covar[["GENE", "avg_all"]], on="GENE", how="inner")
    clusters = enriched_clusters(trait)
    spec = covar.set_index("GENE").loc[genes["GENE"], [f"c{c}" for c in clusters]].to_numpy()
    genes["score"] = (genes["ZSTAT"].to_numpy()[:, None] * spec).max(axis=1)
    genes["rank_pct"] = genes["score"].rank(pct=True)
    return genes, clusters


def matched_permutation(genes, target_symbols, rng):
    """Median rank percentile of the targets against random sets matched on gene size and expression."""
    genes = genes.copy()
    genes["size_bin"] = pd.qcut(genes["NSNPS"], 5, labels=False, duplicates="drop")
    genes["expr_bin"] = pd.qcut(genes["avg_all"], 5, labels=False, duplicates="drop")
    is_target = genes["symbol"].isin(target_symbols)
    observed = genes.loc[is_target, "rank_pct"].median()
    pool = genes.loc[~is_target]
    groups = {k: g["rank_pct"].to_numpy() for k, g in pool.groupby(["size_bin", "expr_bin"])}
    need = genes.loc[is_target].groupby(["size_bin", "expr_bin"]).size().to_dict()
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        draw = [rng.choice(groups[k], size=n, replace=len(groups[k]) < n)
                for k, n in need.items() if k in groups]
        null[i] = np.median(np.concatenate(draw))
    return observed, (np.sum(null >= observed) + 1) / (N_PERM + 1), int(is_target.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True)
    ap.add_argument("--disease", help="benchmark disease label in data/benchmark_targets.tsv")
    args = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)

    genes, clusters = score_genes(args.trait)
    print(f"{args.trait}: {len(genes):,} genes scored from {len(clusters)} clusters")

    if args.disease:
        bt = pd.read_csv("D:/AAN/data/benchmark_targets.tsv", sep="\t")
        bt = bt[(bt["disease"] == args.disease) & (bt["n_targets_of_drug"] <= 3)]
        targets = set(bt["target"])
        obs, p, n_found = matched_permutation(genes, targets, np.random.default_rng(SEED))
        print(f"known targets scored: {n_found} of {len(targets)}; median rank percentile {obs:.3f}; permutation p = {p:.4f}")
        pd.DataFrame([{"disease": args.disease, "n_targets": len(targets), "n_scored": n_found,
                       "median_rank_pct": obs, "perm_p": p, "pass": obs >= 0.80 and p < 0.05}]).to_csv(
            RESULTS / "benchmark.tsv", sep="\t", mode="a", index=False,
            header=not (RESULTS / "benchmark.tsv").exists())
    else:
        top = genes.sort_values("score", ascending=False).head(200)
        top[["GENE", "symbol", "ZSTAT", "P", "score", "rank_pct"]].to_csv(
            RESULTS / f"candidate_genes_{args.trait}.tsv", sep="\t", index=False, float_format="%.5g")
        print(top.head(15)[["symbol", "ZSTAT", "score"]].to_string(index=False))


if __name__ == "__main__":
    main()
