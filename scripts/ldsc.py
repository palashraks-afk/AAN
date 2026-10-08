"""SNP heritability by LD score regression (Bulik-Sullivan et al. 2015), univariate, free intercept.

The reference ldsc package would not build on this Windows machine (it needs pysam), so this is a
small re-implementation of the same regression: chi2_j = 1*a + N_j * h2 * l_j / M, fitted by
iteratively re-weighted least squares with a block jackknife for the standard errors.
It is checked against a simulation in tests/test_ldsc.py.

    python scripts/ldsc.py --name decodeme_gwas_1 --format decodeme --gwas D:/AAN_data/decodeme/gwas_1.regenie.gz
"""
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist
from scipy.stats import norm

LD_DIR = Path("D:/AAN_data/tools/ldsc_ref/eur_w_ld_chr")


def load_ld_scores(ld_dir=LD_DIR):
    frames, m_total = [], 0.0
    for c in range(1, 23):
        frames.append(pd.read_csv(ld_dir / f"{c}.l2.ldscore.gz", sep="\t"))
        m_total += float(open(ld_dir / f"{c}.l2.M_5_50").read().split()[0])
    ld = pd.concat(frames, ignore_index=True)
    return ld[["SNP", "CHR", "BP", "L2"]], m_total


def weights(ld, w_ld, n, m, h2, intercept):
    h2 = min(max(h2, 0.0), 1.0)
    ld = np.maximum(ld, 1.0)
    w_ld = np.maximum(w_ld, 1.0)
    c = h2 * n / m
    het = 1.0 / (2.0 * (intercept + c * ld) ** 2)
    return het / w_ld


def wls(x, y, w):
    sw = np.sqrt(w)
    xw = x * sw[:, None]
    return np.linalg.solve(xw.T @ xw, xw.T @ (y * sw))


def wls_blocks(x, y, w, block_ids, n_blocks):
    """Weighted least squares plus the leave-one-block-out estimates for the jackknife."""
    sw = np.sqrt(w)
    xw = x * sw[:, None]
    yw = y * sw
    xtx = xw.T @ xw
    xty = xw.T @ yw
    est = np.linalg.solve(xtx, xty)
    loo = np.empty((n_blocks, x.shape[1]))
    for b in range(n_blocks):
        rows = block_ids == b
        loo[b] = np.linalg.solve(xtx - xw[rows].T @ xw[rows], xty - xw[rows].T @ yw[rows])
    return est, loo


def jackknife_se(est, loo):
    n = len(loo)
    pseudo = n * est - (n - 1) * loo
    return np.sqrt(pseudo.var(axis=0, ddof=1) / n)


def ld_score_regression(chisq, ld, w_ld, n, m, n_blocks=200, n_iter=2):
    """Return dict with h2, h2_se, intercept, intercept_se, mean_chi2, lambda_gc, n_snps."""
    keep = chisq < max(80.0, 0.001 * n.max())
    chisq, ld, w_ld, n = chisq[keep], ld[keep], w_ld[keep], n[keep]
    nbar = n.mean()
    x = np.column_stack([n * ld / nbar, np.ones(len(chisq))])
    h2 = max((chisq.mean() - 1.0) / np.mean(n * ld / m), 0.0)
    intercept = 1.0
    for _ in range(n_iter):
        w = weights(ld, w_ld, n, m, h2, intercept)
        w = w / w.sum()
        coef = wls(x, chisq, w)
        h2 = m * coef[0] / nbar
        intercept = coef[1]
    w = weights(ld, w_ld, n, m, h2, intercept)
    w = w / w.sum()
    # contiguous blocks along the genome-ordered SNPs
    block_ids = np.minimum((np.arange(len(chisq)) * n_blocks) // len(chisq), n_blocks - 1)
    coef, loo = wls_blocks(x, chisq, w, block_ids, n_blocks)
    se = jackknife_se(coef, loo)
    return {
        "h2_obs": m * coef[0] / nbar, "h2_se": m * se[0] / nbar,
        "intercept": coef[1], "intercept_se": se[1],
        "mean_chi2": float(chisq.mean()), "lambda_gc": float(np.median(chisq) / 0.4549),
        "n_snps": int(len(chisq)),
    }


def liability_scale(h2_obs, h2_se, cases, controls, prevalence):
    """Convert observed-scale h2 to the liability scale (Lee et al. 2011)."""
    p = cases / (cases + controls)
    z = norm.pdf(norm.ppf(1 - prevalence))
    factor = prevalence * (1 - prevalence) / z ** 2 * prevalence * (1 - prevalence) / (p * (1 - p))
    return h2_obs * factor, h2_se * factor


def chisq_from_p(p):
    return chi2_dist.isf(np.clip(p, 1e-300, 1.0), 1)


def prepare(df, ld):
    """Join a GWAS table (rsid, chisq, n) to the LD scores, ordered by genome position."""
    merged = df.merge(ld, left_on="rsid", right_on="SNP", how="inner")
    merged = merged.sort_values(["CHR", "BP"]).reset_index(drop=True)
    return merged


def main():
    import gwas_readers
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--format", choices=["decodeme", "harmonised", "pgc3", "finngen"], required=True)
    ap.add_argument("--n", type=int)
    ap.add_argument("--map", default="D:/AAN_data/derived/decodeme_rsid_map.parquet")
    ap.add_argument("--out", default="D:/AAN_data/derived/ldsc_results.tsv")
    args = ap.parse_args()

    ld, m = load_ld_scores()
    if args.format == "decodeme":
        raw = gwas_readers.read_decodeme(args.gwas, pd.read_parquet(args.map), min_maf=0.0)
    elif args.format == "harmonised":
        raw = gwas_readers.read_harmonised(args.gwas, args.n)
    elif args.format == "finngen":
        raw = gwas_readers.read_finngen(args.gwas, args.n)
    else:
        raw = gwas_readers.read_pgc3(args.gwas)
    raw = raw.drop_duplicates("rsid")
    raw["chisq"] = chisq_from_p(raw["p"].to_numpy())
    merged = prepare(raw, ld)
    print(f"{len(merged):,} SNPs overlap the LD score file (M = {m:,.0f})")
    res = ld_score_regression(merged["chisq"].to_numpy(), merged["L2"].to_numpy(),
                              merged["L2"].to_numpy(), merged["n"].to_numpy().astype(float), m)
    res["h2_z"] = res["h2_obs"] / res["h2_se"]
    res["name"] = args.name
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()})
    row = pd.DataFrame([res])
    out = Path(args.out)
    row.to_csv(out, sep="\t", mode="a", header=not out.exists(), index=False)


if __name__ == "__main__":
    main()
