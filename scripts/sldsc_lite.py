"""C1: partitioned-heritability (S-LDSC style) test of every atlas cluster, written from scratch.

Why: MAGMA and S-LDSC are different statistical models. If a cluster is enriched under both, the result does not
depend on one method's assumptions. This is a simplified S-LDSC: LD scores are computed here from the 1000G EUR
genotypes over HapMap3 SNPs only (the standard pipeline uses all reference SNPs and a 75-annotation baseline), so the
numbers are for testing and ranking clusters, not for quoting enrichment fold-changes.

Model per cluster c, regression SNPs = HapMap3:
    chi2_j = N_j * (t_all * l_all_j + t_expr * l_expr_j + t_c * l_c_j) + a + noise
with l_x the LD score of annotation x (all SNPs; SNPs near any expressed gene; SNPs near the cluster's top-10% genes).
Test: t_c > 0, 200-block jackknife.

    python scripts/sldsc_lite.py --prepare                 # once: LD scores for every annotation
    python scripts/sldsc_lite.py --name decodeme_gwas_1 --format decodeme --gwas <file>
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

import gwas_readers
from ldsc import chisq_from_p, load_ld_scores

REF = Path("C:/AAN_ref/by_chr")
DERIVED = Path("D:/AAN_data/derived")
CACHE = DERIVED / "sldsc"
GENELOC37 = "D:/AAN_data/tools/ref/NCBI37.3.gene.loc"
WINDOW_BP = 1_000_000          # LD window, roughly 1 cM
GENE_FLANK = 100_000
N_BLOCKS = 200
GENO_LUT = np.array([[2, np.nan, 1, 0][(b >> (2 * i)) & 3] for b in range(256) for i in range(4)],
                    dtype=np.float32).reshape(256, 4)


def read_bed_rows(prefix, rows, n_ind):
    row_bytes = (n_ind + 3) // 4
    mm = np.memmap(f"{prefix}.bed", dtype=np.uint8, mode="r", offset=3).reshape(-1, row_bytes)
    packed = mm[rows]
    geno = GENO_LUT[packed].reshape(len(rows), -1)[:, :n_ind]
    return geno


def standardise(geno):
    mean = np.nanmean(geno, axis=1, keepdims=True)
    geno = np.where(np.isnan(geno), mean, geno)
    geno = geno - geno.mean(axis=1, keepdims=True)
    sd = geno.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    return (geno / sd).astype(np.float32)


def merge_windows(starts, ends):
    order = np.argsort(starts)
    s, e = np.asarray(starts)[order], np.asarray(ends)[order]
    out_s, out_e = [s[0]], [e[0]]
    for a, b in zip(s[1:], e[1:]):
        if a <= out_e[-1]:
            out_e[-1] = max(out_e[-1], b)
        else:
            out_s.append(a)
            out_e.append(b)
    return np.array(out_s), np.array(out_e)


def in_windows(pos, starts, ends):
    if len(starts) == 0:
        return np.zeros(len(pos), dtype=bool)
    s, e = merge_windows(starts, ends)
    idx = np.searchsorted(s, pos, side="right") - 1
    ok = idx >= 0
    out = np.zeros(len(pos), dtype=bool)
    out[ok] = pos[ok] <= e[idx[ok]]
    return out


def prepare():
    CACHE.mkdir(parents=True, exist_ok=True)
    ld, _ = load_ld_scores()
    ld = ld.rename(columns={"SNP": "rsid"})
    loc = pd.read_csv(GENELOC37, sep="\t", header=None, names=["GENE", "chr", "start", "end", "strand", "symbol"])
    loc["chr"] = loc["chr"].astype(str)
    sets = {}
    for line in open(DERIVED / "top10_sets.txt"):
        name, genes = line.rstrip("\n").split("\t")
        sets[name] = set(int(g) for g in genes.split())
    names = list(sets)
    expressed = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t", usecols=["GENE"])["GENE"].to_numpy()
    n_ind = sum(1 for _ in open(f"{REF}/g1000_eur_chr1.fam"))

    all_snps, l_all, l_expr, l_clu = [], [], [], []
    for chrom in range(1, 23):
        bim = pd.read_csv(f"{REF}/g1000_eur_chr{chrom}.bim", sep="\t", header=None,
                          names=["chr", "rsid", "cm", "bp", "a1", "a2"])
        bim["row"] = np.arange(len(bim))
        snps = ld[ld["CHR"] == chrom].merge(bim[["rsid", "row"]], on="rsid").sort_values("BP").reset_index(drop=True)
        pos = snps["BP"].to_numpy()
        g = loc[loc["chr"] == str(chrom)].set_index("GENE")
        starts_all = g["start"].to_numpy() - GENE_FLANK
        ends_all = g["end"].to_numpy() + GENE_FLANK
        keep_expr = g.index.isin(expressed)
        A = np.zeros((len(snps), len(names) + 1), dtype=np.float32)
        A[:, 0] = in_windows(pos, starts_all[keep_expr], ends_all[keep_expr])           # expressed genes
        for j, nm in enumerate(names):
            members = g.index.isin(list(sets[nm]))
            A[:, j + 1] = in_windows(pos, starts_all[members], ends_all[members])
        Z = standardise(read_bed_rows(f"{REF}/g1000_eur_chr{chrom}", snps["row"].to_numpy(), n_ind))
        n = Z.shape[1]
        l1 = np.zeros(len(snps), dtype=np.float32)
        lA = np.zeros((len(snps), A.shape[1]), dtype=np.float32)
        B = 1000
        for b0 in range(0, len(snps), B):
            b1 = min(b0 + B, len(snps))
            lo = np.searchsorted(pos, pos[b0] - WINDOW_BP, side="left")
            hi = np.searchsorted(pos, pos[b1 - 1] + WINDOW_BP, side="right")
            r = (Z[b0:b1] @ Z[lo:hi].T) / n
            r2 = r * r - (1 - r * r) / (n - 2)                                  # unbiased r^2
            mask = np.abs(pos[b0:b1, None] - pos[None, lo:hi]) <= WINDOW_BP
            r2 *= mask
            l1[b0:b1] = r2.sum(axis=1)
            lA[b0:b1] = r2 @ A[lo:hi]
        all_snps.append(snps[["rsid", "CHR", "BP"]])
        l_all.append(l1)
        l_expr.append(lA[:, 0])
        l_clu.append(lA[:, 1:])
        print(f"chr{chrom}: {len(snps):,} HapMap3 SNPs", flush=True)
    pd.concat(all_snps, ignore_index=True).to_parquet(CACHE / "snps.parquet")
    np.save(CACHE / "l_all.npy", np.concatenate(l_all))
    np.save(CACHE / "l_expr.npy", np.concatenate(l_expr))
    np.save(CACHE / "l_clusters.npy", np.concatenate(l_clu))
    pd.Series(names).to_csv(CACHE / "cluster_names.txt", index=False, header=False)


def block_wls(X, y, w, n_blocks=N_BLOCKS):
    """Weighted least squares with a block jackknife; returns coefficients and standard errors."""
    sw = np.sqrt(w)[:, None]
    Xw, yw = X * sw, y * sw[:, 0]
    edges = np.linspace(0, len(y), n_blocks + 1).astype(int)       # SNPs are in genome order, blocks are contiguous
    p = X.shape[1]
    xtx = np.zeros((n_blocks, p, p))
    xty = np.zeros((n_blocks, p))
    for b in range(n_blocks):
        sl = slice(edges[b], edges[b + 1])
        xtx[b] = Xw[sl].T @ Xw[sl]
        xty[b] = Xw[sl].T @ yw[sl]
    tot_xx, tot_xy = xtx.sum(axis=0), xty.sum(axis=0)
    est = np.linalg.solve(tot_xx, tot_xy)
    loo = np.array([np.linalg.solve(tot_xx - xtx[b], tot_xy - xty[b]) for b in range(n_blocks)])
    pseudo = n_blocks * est - (n_blocks - 1) * loo
    return est, np.sqrt(pseudo.var(axis=0, ddof=1) / n_blocks)


def run_trait(name, raw):
    snps = pd.read_parquet(CACHE / "snps.parquet")
    l_all = np.load(CACHE / "l_all.npy")
    l_expr = np.load(CACHE / "l_expr.npy")
    l_clu = np.load(CACHE / "l_clusters.npy", mmap_mode="r")
    raw = raw.drop_duplicates("rsid").copy()
    raw["chisq"] = chisq_from_p(raw["p"].to_numpy())
    idx = snps.reset_index().merge(raw[["rsid", "chisq", "n"]], on="rsid")
    idx = idx[idx["chisq"] < max(80.0, 0.001 * idx["n"].max())]
    rows = idx["index"].to_numpy()
    y = idx["chisq"].to_numpy()
    n = idx["n"].to_numpy().astype(float)
    la, le = l_all[rows], l_expr[rows]
    print(f"{name}: {len(rows):,} HapMap3 SNPs in the regression", flush=True)

    names = pd.read_csv(CACHE / "cluster_names.txt", header=None)[0].tolist()
    out = []
    base = np.column_stack([n * la, n * le, np.ones(len(y))])
    chunk_start, chunk = -1, None
    for j, nm in enumerate(names):
        if j >= chunk_start + 50 or chunk is None:                      # read the cluster LD scores 50 columns at a time
            chunk_start = j
            chunk = np.asarray(l_clu[:, j:j + 50])[rows]
        X = np.column_stack([base[:, 0], base[:, 1], n * chunk[:, j - chunk_start], base[:, 2]])
        w = 1.0 / np.maximum(la, 1.0)
        for _ in range(2):                                      # reweight by the fitted mean (var(chi2) ~ 2 mean^2)
            est, _se = block_wls(X, y, w / w.sum(), n_blocks=20)
            fitted = np.maximum(X @ est, 1.0)
            w = 1.0 / (np.maximum(la, 1.0) * fitted ** 2)
        est, se = block_wls(X, y, w / w.sum())
        z = est[2] / se[2]
        out.append({"cluster": nm, "tau": est[2], "se": se[2], "z": z, "p_one_sided": norm.sf(z)})
    res = pd.DataFrame(out)
    res.to_csv(Path("D:/AAN/results") / f"sldsc_{name}.tsv", sep="\t", index=False, float_format="%.5g")
    print(res.sort_values("p_one_sided").head(8).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--name")
    ap.add_argument("--format", choices=["decodeme", "harmonised", "pgc3", "finngen"])
    ap.add_argument("--gwas")
    ap.add_argument("--n", type=int)
    a = ap.parse_args()
    if a.prepare:
        prepare()
    else:
        if a.format == "decodeme":
            raw = gwas_readers.read_decodeme(a.gwas, pd.read_parquet(DERIVED / "decodeme_rsid_map.parquet"), min_maf=0.0)
        elif a.format == "harmonised":
            raw = gwas_readers.read_harmonised(a.gwas, a.n)
        elif a.format == "finngen":
            raw = gwas_readers.read_finngen(a.gwas, a.n)
        else:
            raw = gwas_readers.read_pgc3(a.gwas)
        run_trait(a.name, raw)
