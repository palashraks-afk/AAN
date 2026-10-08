"""Does higher expression of a candidate gene go with higher or lower ME/CFS risk? (SMR without HEIDI)

For each candidate gene and each GTEx v10 brain tissue where it is an eGene (q < 0.05), take the lead cis-eQTL
variant, look it up in the ME/CFS GWAS, and compute the summary-data-based Mendelian randomisation statistic
(Zhu et al. 2016): T = z_gwas^2 * z_eqtl^2 / (z_gwas^2 + z_eqtl^2), chi-square with 1 df. The sign of
beta_gwas / beta_eqtl tells the direction: positive means more expression, more risk.

Limits that must be stated with any result: no HEIDI test, so a shared signal cannot be told apart from two nearby
causal variants; one instrument per gene and tissue; GTEx donors are mostly European but not DecodeME.

    python scripts/smr_lite.py --candidates results/candidate_genes_decodeme_gwas_1.tsv --top 200
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

GTEX = Path("D:/AAN_data/translational/gtex_brain/GTEx_Analysis_v10_eQTL_updated")
GWAS = "D:/AAN_data/decodeme/gwas_1.regenie.gz"
QC_LIST = "D:/AAN_data/decodeme/gwas_qced.var.gz"
RESULTS = Path("D:/AAN/results")


def load_egenes(symbols):
    frames = []
    for path in sorted(glob.glob(str(GTEX / "*.v10.eGenes.txt.gz"))):
        tissue = Path(path).name.replace(".v10.eGenes.txt.gz", "")
        e = pd.read_csv(path, sep="\t", usecols=["gene_name", "chr", "variant_pos", "ref", "alt", "af",
                                                  "slope", "slope_se", "qval", "variant_id"])
        e = e[e["gene_name"].isin(symbols) & (e["qval"] < 0.05)].copy()
        e["tissue"] = tissue
        frames.append(e)
    out = pd.concat(frames, ignore_index=True)
    out["chrom"] = out["chr"].str.replace("chr", "", regex=False)
    return out[out["chrom"].isin([str(i) for i in range(1, 23)])]


def lookup_gwas(eqtl):
    needed = set(zip(eqtl["chrom"], eqtl["variant_pos"]))
    passed = set(pd.read_csv(QC_LIST, header=None)[0])
    hits = []
    for chunk in pd.read_csv(GWAS, sep=r"\s+", usecols=["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "BETA", "SE"],
                             dtype={"CHROM": str}, chunksize=1_000_000):
        key = list(zip(chunk["CHROM"], chunk["GENPOS"]))
        hit = chunk[[k in needed for k in key]]
        hits.append(hit[hit["ID"].isin(passed)])
    return pd.concat(hits, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--top", type=int, default=200)
    args = ap.parse_args()

    cand = pd.read_csv(args.candidates, sep="\t").head(args.top)
    eqtl = load_egenes(set(cand["symbol"]))
    print(f"{len(eqtl)} gene x tissue lead eQTLs for {eqtl['gene_name'].nunique()} of {len(cand)} candidates")
    gw = lookup_gwas(eqtl)

    m = eqtl.merge(gw, left_on=["chrom", "variant_pos"], right_on=["CHROM", "GENPOS"], how="inner")
    same_pair = [{r, a} == {a0, a1} for r, a, a0, a1 in zip(m["ref"], m["alt"], m["ALLELE0"], m["ALLELE1"])]
    m = m[same_pair].copy()
    # eQTL slope is for the alt allele, the GWAS beta is for ALLELE1
    sign = np.where(m["alt"] == m["ALLELE1"], 1.0, -1.0)
    m["beta_gwas"] = m["BETA"] * sign
    m["z_eqtl"] = m["slope"] / m["slope_se"]
    m["z_gwas"] = m["beta_gwas"] / m["SE"]
    t = (m["z_gwas"] ** 2 * m["z_eqtl"] ** 2) / (m["z_gwas"] ** 2 + m["z_eqtl"] ** 2)
    m["smr_chi2"] = t
    m["smr_p"] = chi2.sf(t, 1)
    m["direction"] = np.where(m["beta_gwas"] / m["slope"] > 0, "higher expression, higher risk",
                              "higher expression, lower risk")
    m["strong_instrument"] = m["z_eqtl"] ** 2 > 28.7        # p < 5e-8, the SMR convention
    keep = ["gene_name", "tissue", "variant_id", "af", "slope", "slope_se", "beta_gwas", "SE", "z_eqtl", "z_gwas",
            "smr_chi2", "smr_p", "direction", "strong_instrument"]
    out = m[keep].rename(columns={"SE": "se_gwas"}).sort_values("smr_p")
    out.to_csv(RESULTS / "smr_lite_decodeme_gwas_1.tsv", sep="\t", index=False, float_format="%.5g")
    n_genes = out["gene_name"].nunique()
    print(f"{n_genes} genes tested; Bonferroni threshold p < {0.05 / max(n_genes, 1):.2e}")
    print(out.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
