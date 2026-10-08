"""Run MAGMA (SNP -> gene -> cell-type gene-property) for one GWAS.

    python scripts/run_magma.py --name decodeme_gwas_1 --format decodeme --gwas D:/AAN_data/decodeme/gwas_1.regenie.gz
    python scripts/run_magma.py --name alzheimer --format harmonised --gwas D:/AAN_data/panel/alzheimer.h.tsv.gz --n 487511

Steps: write the p-value and SNP-location files, annotate SNPs to genes (window 35 kb up, 10 kb
down), gene analysis against the 1000G EUR panel, then the gene-property test for every cluster,
conditioning on the two expression covariates.
"""
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

import gwas_readers

TOOLS = Path("D:/AAN_data/tools")
MAGMA = TOOLS / "magma" / "magma.exe"
REF = Path("C:/AAN_ref/g1000_eur")   # copied to the SSD, the D: drive is a slow hard disk
GENELOC = {"GRCh38": TOOLS / "ref" / "NCBI38.gene.loc", "GRCh37": TOOLS / "ref" / "NCBI37.3.gene.loc"}
DERIVED = Path("D:/AAN_data/derived")


def run(args):
    print(" ".join(str(a) for a in args), flush=True)
    subprocess.run([str(a) for a in args], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--format", choices=["decodeme", "harmonised", "pgc3"], required=True)
    ap.add_argument("--n", type=int, help="total sample size, used when the file has no per-SNP N")
    ap.add_argument("--map", default=str(DERIVED / "decodeme_rsid_map.parquet"))
    ap.add_argument("--covar", default=str(DERIVED / "gene_covar.txt"))
    ap.add_argument("--condition", default="avg_all,avg_neuron")
    ap.add_argument("--jobs", type=int, default=6, help="chromosomes analysed in parallel")
    args = ap.parse_args()

    build = "GRCh37" if args.format == "pgc3" else "GRCh38"
    workdir = DERIVED / "magma" / args.name
    workdir.mkdir(parents=True, exist_ok=True)

    if args.format == "decodeme":
        g = gwas_readers.read_decodeme(args.gwas, pd.read_parquet(args.map))
    elif args.format == "harmonised":
        g = gwas_readers.read_harmonised(args.gwas, args.n)
    else:
        g = gwas_readers.read_pgc3(args.gwas)
    g = g.drop_duplicates("rsid")
    print(f"{len(g):,} SNPs read ({build})", flush=True)

    g[["rsid", "p", "n"]].to_csv(workdir / "pval.txt", sep="\t", index=False,
                                 header=["SNP", "P", "N"], float_format="%.6g")
    g[["rsid", "chrom", "pos"]].to_csv(workdir / "snploc.txt", sep="\t", index=False, header=False)

    run([MAGMA, "--annotate", "window=35,10", "--snp-loc", workdir / "snploc.txt",
         "--gene-loc", GENELOC[build], "--out", workdir / "annot"])
    # one job per chromosome, a few at a time: much faster than a single run and it keeps memory down
    def gene_analysis(chrom):
        run([MAGMA, "--bfile", REF, "--pval", workdir / "pval.txt", "ncol=N",
             "--gene-annot", workdir / "annot.genes.annot", "--batch", chrom, "chr",
             "--out", workdir / "genes"])

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(gene_analysis, [str(c) for c in range(1, 23)]))
    run([MAGMA, "--merge", workdir / "genes", "--out", workdir / "genes"])

    run([MAGMA, "--gene-results", workdir / "genes.genes.raw", "--gene-covar", args.covar,
         "--model", f"condition-hide={args.condition}", "direction=pos",
         "--out", workdir / "cluster_prop"])


if __name__ == "__main__":
    main()
