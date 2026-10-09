"""L1 (PREREG v1.5): how many independent cell-type signals does a trait have? Forward selection with MAGMA conditioning.

    python scripts/independent_signals.py --name schizophrenia

Step k conditions on mean expression plus the clusters already chosen, retests every other cluster, and adds the best one
if its FDR is below 0.05, up to 5 clusters. Run it on the controls first (schizophrenia, Alzheimer's, rheumatoid arthritis).
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

import calibrate

DERIVED = Path("D:/AAN_data/derived")
RESULTS = Path("D:/AAN/results")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")
MAX_SIGNALS = 5


def run_step(workdir, selected, tag):
    cond = ",".join(["avg_all"] + [f"c{c}" for c in selected])
    out = workdir / f"indep_{tag}"
    subprocess.run([str(MAGMA), "--gene-results", str(workdir / "genes.genes.raw"), "--gene-covar",
                    str(DERIVED / "gene_covar.txt"), "--model", f"condition-hide={cond}", "direction=pos",
                    "--out", str(out)], check=True, stdout=subprocess.DEVNULL)
    df = pd.read_csv(f"{out}.gsa.out.txt", comment="#", sep=r"\s+")
    df = df[df["VARIABLE"].str.match(r"^c\d+$")].copy()
    df["cluster_id"] = df["VARIABLE"].str[1:].astype(int)
    return df.set_index("cluster_id")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    args = ap.parse_args()
    workdir = DERIVED / "magma" / args.name
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    unconditional = calibrate.load_model(args.name, "A")

    selected, rows = [], []
    for step in range(MAX_SIGNALS):
        res = run_step(workdir, selected, f"{step}")
        res = res.drop(index=[c for c in selected if c in res.index])
        res["fdr"] = multipletests(res["P"], method="fdr_bh")[1]
        best = res["P"].idxmin()
        if res.loc[best, "fdr"] >= 0.05:
            print(f"step {step + 1}: nothing else passes FDR < 0.05, stopping")
            break
        rows.append({"step": step + 1, "cluster_id": int(best), "cluster": ann.loc[best, "cluster_name"],
                     "supercluster": ann.loc[best, "supercluster"],
                     "z_conditional": calibrate.z_from_p(res.loc[best, "P"]),
                     "z_unconditional": calibrate.z_from_p(unconditional.loc[best, "P"]),
                     "fdr_conditional": res.loc[best, "fdr"]})
        selected.append(int(best))
        print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / f"independent_signals_{args.name}.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"{args.name}: {len(out)} independent cell-type signal(s)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
