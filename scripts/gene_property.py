"""MAGMA gene-property test of every atlas cluster, for one trait and one pre-registered model.

Models (see preregistration/PREREG_v1.md):
    A  condition on avg_all
    B  condition on avg_all and avg_neuron
    C  condition on avg_all, avg_neuron and the gene-level z of depression, BMI and insomnia

    python scripts/gene_property.py --name decodeme_gwas_1 --model A
"""
import argparse
import subprocess
from pathlib import Path

DERIVED = Path("D:/AAN_data/derived")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")

MODELS = {
    "A": ["avg_all"],
    "B": ["avg_all", "avg_neuron"],
    "C": ["avg_all", "avg_neuron", "z_depression_broad", "z_bmi", "z_insomnia"],
}


def run_model(name, model):
    workdir = DERIVED / "magma" / name
    covar = DERIVED / ("gene_covar_C.txt" if model == "C" else "gene_covar.txt")
    cmd = [MAGMA, "--gene-results", workdir / "genes.genes.raw", "--gene-covar", covar,
           "--model", "condition-hide=" + ",".join(MODELS[model]), "direction=pos",
           "--out", workdir / f"cluster_{model}"]
    print(" ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--model", choices=sorted(MODELS), required=True)
    a = ap.parse_args()
    run_model(a.name, a.model)
