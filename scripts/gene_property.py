"""MAGMA gene-property and gene-set tests on the merged gene results of one trait.

Main cell-type models (see preregistration/PREREG_v1.md):
    A  condition on avg_all
    B  condition on avg_all and avg_neuron
    C  condition on avg_all, avg_neuron and the gene-level z of depression, BMI and insomnia

Discovery layer (preregistration/PREREG_v1.2_discovery_layer.md):
    groups   transmitter-class and neuron-subtype groups, forms A and B
    hpa      Human Protein Atlas cell types from the whole body, form A
    reactome Reactome pathways, conditioning on avg_all and avg_neuron; with --top-cluster also on that cluster

    python scripts/gene_property.py --name decodeme_gwas_1 --model A
    python scripts/gene_property.py --name decodeme_gwas_1 --layer groups
    python scripts/gene_property.py --name decodeme_gwas_1 --layer reactome --top-cluster 253
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd

DERIVED = Path("D:/AAN_data/derived")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")

MODELS = {
    "A": ["avg_all"],
    "B": ["avg_all", "avg_neuron"],
    "C": ["avg_all", "avg_neuron", "z_depression_broad", "z_bmi", "z_insomnia"],
}


def magma(cmd):
    print(" ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def run_model(name, model):
    workdir = DERIVED / "magma" / name
    covar = DERIVED / ("gene_covar_C.txt" if model == "C" else "gene_covar.txt")
    magma([MAGMA, "--gene-results", workdir / "genes.genes.raw", "--gene-covar", covar,
           "--model", "condition-hide=" + ",".join(MODELS[model]), "direction=pos",
           "--out", workdir / f"cluster_{model}"])


def run_groups(name):
    workdir = DERIVED / "magma" / name
    for form, cond in (("A", ["avg_all"]), ("B", ["avg_all", "avg_neuron"])):
        magma([MAGMA, "--gene-results", workdir / "genes.genes.raw", "--gene-covar",
               DERIVED / "gene_covar_groups.txt", "--model", "condition-hide=" + ",".join(cond), "direction=pos",
               "--out", workdir / f"groups_{form}"])


def run_hpa(name):
    workdir = DERIVED / "magma" / name
    magma([MAGMA, "--gene-results", workdir / "genes.genes.raw", "--gene-covar", DERIVED / "gene_covar_hpa.txt",
           "--model", "condition-hide=avg_all", "direction=pos", "--out", workdir / "hpa_A"])


def run_regions(name):
    workdir = DERIVED / "magma" / name
    magma([MAGMA, "--gene-results", workdir / "genes.genes.raw", "--gene-covar", DERIVED / "gene_covar_regions.txt",
           "--model", "condition-hide=avg_all", "direction=pos", "--out", workdir / "regions_A"])


def run_reactome(name, top_cluster=None):
    workdir = DERIVED / "magma" / name
    base = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t", usecols=["GENE", "avg_all", "avg_neuron"]
                       + ([f"c{top_cluster}"] if top_cluster is not None else []))
    cond = ["avg_all", "avg_neuron"]
    out_prefix = "reactome_A"
    if top_cluster is not None:
        base = base.rename(columns={f"c{top_cluster}": "top_cluster_spec"})
        cond.append("top_cluster_spec")
        out_prefix = "reactome_C"
    covar = workdir / f"reactome_covar_{out_prefix}.txt"
    base.to_csv(covar, sep="\t", index=False, float_format="%.5f")
    magma([MAGMA, "--gene-results", workdir / "genes.genes.raw", "--set-annot", DERIVED / "reactome.sets",
           "--gene-covar", covar, "--model", "condition-hide=" + ",".join(cond), "analyse=sets", "direction=pos",
           "--out", workdir / out_prefix])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--model", choices=sorted(MODELS))
    ap.add_argument("--layer", choices=["groups", "hpa", "regions", "reactome"])
    ap.add_argument("--top-cluster", type=int)
    a = ap.parse_args()
    if a.model:
        run_model(a.name, a.model)
    if a.layer == "groups":
        run_groups(a.name)
    elif a.layer == "hpa":
        run_hpa(a.name)
    elif a.layer == "regions":
        run_regions(a.name)
    elif a.layer == "reactome":
        run_reactome(a.name, a.top_cluster)
