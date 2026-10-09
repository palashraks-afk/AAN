"""C2 and C3 of preregistration/PREREG_v1.3_confirmation.md for one trait.

    python scripts/confirm_robustness.py --name decodeme_gwas_1 --step locus      # locus drop + leave-one-chromosome-out
    python scripts/confirm_robustness.py --name decodeme_gwas_1 --step permute    # gene-level permutation null

Every step reruns MAGMA's gene-property model A on a modified covariate file. MAGMA drops genes that are missing from
the covariate file, which is how genes are removed from a run.
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import calibrate

DERIVED = Path("D:/AAN_data/derived")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")
GENELOC = "D:/AAN_data/tools/ref/NCBI38.gene.loc"
N_PERM_TOTAL = 1000          # preregistered 5,000; reduced to 1,000 for runtime, logged in DEVIATIONS.md
PERM_BATCH = 100             # permutations per MAGMA run (x 20 clusters = 2,000 columns)
TOP_CLUSTERS = 20
SEED = 20261008


def gene_positions():
    loc = pd.read_csv(GENELOC, sep="\t", header=None, names=["GENE", "chr", "start", "end", "strand", "symbol"])
    loc["chr"] = loc["chr"].astype(str)
    return loc


def run_model_a(workdir, covar_path, out_prefix, condition="avg_all"):
    cmd = [MAGMA, "--gene-results", workdir / "genes.genes.raw", "--gene-covar", covar_path,
           "--model", f"condition-hide={condition}", "direction=pos", "--out", workdir / out_prefix]
    subprocess.run([str(c) for c in cmd], check=True, stdout=subprocess.DEVNULL)


def locus_genes(workdir, window=1_000_000):
    """Genes within 1 Mb of any SNP with p < 5e-8."""
    p = pd.read_csv(workdir / "pval.txt", sep="\t")
    loc = pd.read_csv(workdir / "snploc.txt", sep="\t", header=None, names=["SNP", "chr", "pos"])
    hits = p.merge(loc, on="SNP")
    hits = hits[hits["P"] < 5e-8]
    genes = gene_positions()
    drop = set()
    for chrom, h in hits.groupby(hits["chr"].astype(str)):
        pos = h["pos"].to_numpy()
        g = genes[genes["chr"] == chrom]
        for gid, s, e in zip(g["GENE"], g["start"], g["end"]):
            if ((pos > s - window) & (pos < e + window)).any():
                drop.add(gid)
    return drop, len(hits)


def step_locus(name):
    workdir = DERIVED / "magma" / name
    covar = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t")
    drop, n_hits = locus_genes(workdir)
    print(f"{n_hits} genome-wide significant SNPs; dropping {len(drop & set(covar['GENE']))} genes within 1 Mb")
    path = workdir / "covar_locusdrop.txt"
    covar[~covar["GENE"].isin(drop)].to_csv(path, sep="\t", index=False, float_format="%.5f")
    run_model_a(workdir, path, "cluster_locusdrop")

    genes = gene_positions().set_index("GENE")["chr"]
    covar["chr"] = covar["GENE"].map(genes)
    for c in range(1, 23):
        sub = covar[covar["chr"] != str(c)].drop(columns="chr")
        path = workdir / f"covar_loco_{c}.txt"
        sub.to_csv(path, sep="\t", index=False, float_format="%.5f")
        run_model_a(workdir, path, f"cluster_loco_chr{c}")
        path.unlink()
        print("left out chromosome", c, flush=True)


def beta_std(workdir, prefix):
    base = workdir / f"{prefix}.gsa.out.txt"
    df = pd.read_csv(base, comment="#", sep=r"\s+")
    return df[df["VARIABLE"].str.match(r"^c\d+$")].set_index("VARIABLE")["BETA_STD"]


def step_permute(name):
    rng = np.random.default_rng(SEED)
    workdir = DERIVED / "magma" / name
    covar = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t")
    top = calibrate.load_model(name, "A")["P"].nsmallest(TOP_CLUSTERS).index.tolist()
    cols = [f"c{c}" for c in top]
    observed = beta_std(workdir, "cluster_A")[cols]

    # genes are shuffled only among genes of similar size and expression
    genes = pd.read_csv(workdir / "genes.genes.out", sep="\t").set_index("GENE")["NSNPS"]
    covar["nsnps"] = covar["GENE"].map(genes)
    covar = covar.dropna(subset=["nsnps"]).reset_index(drop=True)
    covar["bin"] = (pd.qcut(covar["nsnps"], 5, labels=False, duplicates="drop") * 5
                    + pd.qcut(covar["avg_all"], 5, labels=False, duplicates="drop"))
    groups = [np.where(covar["bin"].to_numpy() == b)[0] for b in sorted(covar["bin"].unique())]
    spec = covar[cols].to_numpy()

    exceed = pd.Series(0, index=cols)
    done = 0
    while done < N_PERM_TOTAL:
        block = {"GENE": covar["GENE"].to_numpy(), "avg_all": covar["avg_all"].to_numpy()}
        for k in range(PERM_BATCH):
            perm = np.arange(len(covar))
            for idx in groups:
                perm[idx] = rng.permutation(idx)
            for j, c in enumerate(cols):
                block[f"{c}_p{k}"] = spec[perm, j]
        path = workdir / "covar_perm.txt"
        pd.DataFrame(block).to_csv(path, sep="\t", index=False, float_format="%.4f")
        run_model_a(workdir, path, "perm_batch")
        b = beta_std(workdir, "perm_batch")
        for c in cols:
            vals = b[[f"{c}_p{k}" for k in range(PERM_BATCH)]].to_numpy()
            exceed[c] += int((vals >= observed[c]).sum())
        done += PERM_BATCH
        print(f"{done}/{N_PERM_TOTAL} permutations", flush=True)
        path.unlink()
    emp = (exceed + 1) / (N_PERM_TOTAL + 1)
    out = pd.DataFrame({"cluster": top, "observed_beta_std": observed.values, "empirical_p": emp.values})
    out.to_csv(Path("D:/AAN/results") / f"permutation_{name}.tsv", sep="\t", index=False, float_format="%.5g")
    print(out.to_string(index=False))


def significant_loci(workdir, window=1_000_000):
    """Merge SNPs with p < 5e-8 into loci (SNPs closer than 1 Mb join the same locus)."""
    p = pd.read_csv(workdir / "pval.txt", sep="	")
    loc = pd.read_csv(workdir / "snploc.txt", sep="	", header=None, names=["SNP", "chr", "pos"])
    hits = p.merge(loc, on="SNP")
    hits = hits[hits["P"] < 5e-8].copy()
    hits["chr"] = hits["chr"].astype(str)
    loci = []
    for chrom, h in hits.groupby("chr"):
        pos = np.sort(h["pos"].to_numpy())
        start = prev = pos[0]
        for x in pos[1:]:
            if x - prev > window:
                loci.append((chrom, start, prev))
                start = x
            prev = x
        loci.append((chrom, start, prev))
    return loci


def step_perlocus(name):
    workdir = DERIVED / "magma" / name
    covar = pd.read_csv(DERIVED / "gene_covar.txt", sep="	")
    genes = gene_positions()
    base = calibrate.load_model(name, "A")
    top = base["P"].nsmallest(3).index.tolist()
    rows = [{"locus": "none (full data)", **{f"z_c{c}": float(calibrate.z_from_p(base.loc[c, "P"])) for c in top}}]
    for chrom, start, end in significant_loci(workdir):
        g = genes[genes["chr"] == chrom]
        drop = set(g.loc[(g["end"] > start - 1_000_000) & (g["start"] < end + 1_000_000), "GENE"])
        path = workdir / "covar_perlocus.txt"
        covar[~covar["GENE"].isin(drop)].to_csv(path, sep="	", index=False, float_format="%.5f")
        run_model_a(workdir, path, "perlocus_tmp")
        m = pd.read_csv(workdir / "perlocus_tmp.gsa.out.txt", comment="#", sep=r"\s+").set_index("VARIABLE")
        rows.append({"locus": f"chr{chrom}:{start}-{end} ({len(drop)} genes dropped)",
                     **{f"z_c{c}": float(calibrate.z_from_p(m.loc[f"c{c}", "P"])) for c in top}})
        print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(Path("D:/AAN/results") / f"perlocus_{name}.tsv", sep="	", index=False, float_format="%.4g")
    print(out.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--step", choices=["locus", "permute", "perlocus"], required=True)
    a = ap.parse_args()
    {"locus": step_locus, "permute": step_permute, "perlocus": step_perlocus}[a.step](a.name)
