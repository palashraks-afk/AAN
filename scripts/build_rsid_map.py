"""Map DecodeME variants (GRCh38, ids like 5:90000126:G:A) to rsIDs.

The 1000G reference that MAGMA and LDSC use is GRCh37 with rsIDs, so I lift each DecodeME
position back to GRCh37 and match on chromosome, position and the pair of alleles.
Only plain SNPs are kept. Output is one parquet table that every later step reuses.
"""
import argparse
from pathlib import Path

import pandas as pd
from pyliftover import LiftOver

AUTOSOMES = {str(i) for i in range(1, 23)}


def load_variants(gwas_path):
    cols = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1"]
    df = pd.read_csv(gwas_path, sep=r"\s+", usecols=cols,
                     dtype={"CHROM": str, "GENPOS": int}, engine="c")
    df = df[df["CHROM"].isin(AUTOSOMES)]
    df = df[(df["ALLELE0"].str.len() == 1) & (df["ALLELE1"].str.len() == 1)]
    return df.reset_index(drop=True)


def lift_positions(df):
    lo = LiftOver("hg38", "hg19")
    pairs = df[["CHROM", "GENPOS"]].drop_duplicates()
    out = {}
    for chrom, pos in zip(pairs["CHROM"], pairs["GENPOS"]):
        hit = lo.convert_coordinate("chr" + chrom, pos - 1)  # liftover is 0-based
        if hit:
            new_chrom, new_pos = hit[0][0], hit[0][1] + 1
            if new_chrom == "chr" + chrom:
                out[(chrom, pos)] = new_pos
    keys = list(zip(df["CHROM"], df["GENPOS"]))
    return pd.Series([out.get(k, -1) for k in keys], index=df.index)


def load_reference(bim_path):
    bim = pd.read_csv(bim_path, sep="\t", header=None,
                      names=["CHROM", "rsid", "cm", "POS37", "A1", "A2"],
                      dtype={"CHROM": str}, usecols=[0, 1, 3, 4, 5])
    bim = bim[bim["CHROM"].isin(AUTOSOMES)]
    bim = bim[(bim["A1"].str.len() == 1) & (bim["A2"].str.len() == 1)]
    bim["allele_pair"] = [a + b if a < b else b + a for a, b in zip(bim["A1"], bim["A2"])]
    return bim[["CHROM", "POS37", "rsid", "allele_pair"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", default="D:/AAN_data/decodeme/gwas_1.regenie.gz")
    ap.add_argument("--bim", default="D:/AAN_data/tools/ref/g1000_eur.bim")
    ap.add_argument("--out", default="D:/AAN_data/derived/decodeme_rsid_map.parquet")
    args = ap.parse_args()

    variants = load_variants(args.gwas)
    print("snps in gwas:", len(variants))

    variants["POS37"] = lift_positions(variants)
    variants = variants[variants["POS37"] > 0]
    print("lifted to GRCh37:", len(variants))

    variants["allele_pair"] = [a + b if a < b else b + a
                               for a, b in zip(variants["ALLELE0"], variants["ALLELE1"])]

    ref = load_reference(args.bim)
    print("reference snps:", len(ref))

    merged = variants.merge(ref, on=["CHROM", "POS37", "allele_pair"], how="inner")
    # a position can match more than one rsID, keep the first so each variant appears once
    merged = merged.drop_duplicates(subset="ID", keep="first")
    print("matched to an rsID:", len(merged), f"({len(merged) / len(variants):.1%})")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    merged[["ID", "rsid", "CHROM", "GENPOS", "POS37", "ALLELE0", "ALLELE1"]].to_parquet(args.out)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
