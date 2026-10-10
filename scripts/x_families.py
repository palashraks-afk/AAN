"""PREREG v1.7 F1: eight pre-fixed HGNC gene families, MAGMA competitive set test conditioned on expression.

    python scripts/x_families.py --names decodeme_gwas_1 decodeme_gwas_2 alzheimer schizophrenia height
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

DERIVED = Path("D:/AAN_data/derived")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")
HGNC = Path("D:/AAN_data/translational/hgnc_complete_set.txt")
RESULTS = Path("D:/AAN/results")
FAMILIES = {  # fixed in PREREG v1.7 before any result
    "ion_channel": "channel", "gpcr": "G protein-coupled receptor", "neuropeptide": "neuropeptide|peptide hormone",
    "solute_carrier": "solute carrier", "kinase": "kinase", "nuclear_receptor": "nuclear receptor",
    "cadherin": "cadherin", "immunoglobulin": "immunoglobulin",
}


def build_sets():
    d = pd.read_csv(HGNC, sep="\t", low_memory=False, usecols=["entrez_id", "gene_group"]).dropna()
    d["entrez_id"] = d["entrez_id"].astype(int)
    path = DERIVED / "families.sets"
    with open(path, "w") as fh:
        for name, pat in FAMILIES.items():
            genes = sorted(set(d.loc[d["gene_group"].str.contains(pat, case=False, regex=True), "entrez_id"]))
            assert len(genes) >= 30, name
            fh.write(name + " " + " ".join(map(str, genes)) + "\n")
    return path


def run(name, sets):
    workdir = DERIVED / "magma" / name
    covar = workdir / "families_covar.txt"
    pd.read_csv(DERIVED / "gene_covar.txt", sep="\t", usecols=["GENE", "avg_all", "avg_neuron"]).to_csv(
        covar, sep="\t", index=False, float_format="%.5f")
    subprocess.run([str(MAGMA), "--gene-results", str(workdir / "genes.genes.raw"), "--set-annot", str(sets),
                    "--gene-covar", str(covar), "--model", "condition-hide=avg_all,avg_neuron", "analyse=sets",
                    "direction=pos", "--out", str(workdir / "families")], check=True)


def read(name):
    for suf in (".gsa.out", ".gsa.out.txt"):
        f = DERIVED / "magma" / name / f"families{suf}"
        if f.exists():
            df = pd.read_csv(f, comment="#", sep=r"\s+")
            df = df[df["TYPE"] == "SET"] if "TYPE" in df.columns else df
            return df.set_index("VARIABLE")[["NGENES", "BETA", "SE", "P"]]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    a = ap.parse_args()
    sets = build_sets()
    rows = []
    for n in a.names:
        run(n, sets)
        t = read(n)
        t["trait"] = n
        t["fdr"] = multipletests(t["P"], method="fdr_bh")[1]
        rows.append(t.reset_index())
    out = pd.concat(rows)
    out.to_csv(RESULTS / "x_F1_families.tsv", sep="\t", index=False, float_format="%.5g")
    print(out.to_string())
