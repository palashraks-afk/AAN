"""PREREG v1.11 R1: FinnGen R13 R18_MALAI_FATIG as an independent fatigue cohort.

    python scripts/x_finngen_fatigue.py
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import calibrate

HERE = Path(__file__).resolve().parent
DERIVED = Path("D:/AAN_data/derived")
RESULTS = Path("D:/AAN/results")
FILE = Path("D:/AAN_data/panel/related/finngen_R13_R18_MALAI_FATIG.gz")
NAME = "finngen_malaise_fatigue"
N_TOTAL = 480_000
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]
ELEVEN = SIX + ["DLIT_136", "DLIT_151", "ULIT_126", "DLIT_148", "DLIT_147"]


def run(cmd):
    print(" ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], cwd=HERE).returncode == 0


def main():
    f = DERIVED / "ldsc_results.tsv"
    have = f.exists() and NAME in pd.read_csv(f, sep="\t")["name"].values
    if not have:
        run([sys.executable, "-u", "ldsc.py", "--name", NAME, "--format", "finngen", "--gwas", FILE, "--n", N_TOTAL])
    t = pd.read_csv(f, sep="\t").drop_duplicates("name", keep="last").set_index("name")
    z = float(t.loc[NAME, "h2_z"])
    print("heritability z:", z)
    if z < 4:
        print("excluded by the heritability gate; no cell-type analysis (pre-registered)")
        return
    if not any((DERIVED / "magma" / NAME / f"cluster_A.gsa.out{s}").exists() for s in (".txt", "")):
        run([sys.executable, "-u", "run_magma.py", "--name", NAME, "--format", "finngen", "--gwas", FILE, "--n", N_TOTAL, "--jobs", 8])
    m = calibrate.load_model(NAME, "A")
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    m["z"] = calibrate.z_from_p(m["P"].to_numpy())
    m = m.join(ann[["cluster_name"]])
    six = m[m["cluster_name"].isin(SIX)]
    obs = six["z"].mean()
    rng = np.random.default_rng(20261010)
    null = np.array([rng.choice(m["z"].to_numpy(), 6, replace=False).mean() for _ in range(10000)])
    p = (1 + (null >= obs).sum()) / 10001
    print(f"mean z of the six specific clusters: {obs:.3f}; random-set 95th percentile {np.percentile(null, 95):.3f}; permutation p = {p:.4f}")
    print("clusters with z > 0:", int((six["z"] > 0).sum()), "of", len(six))
    print(six[["cluster_name", "z", "P"]].to_string())
    eleven = m[m["cluster_name"].isin(ELEVEN)][["cluster_name", "z", "P"]]
    eleven.to_csv(RESULTS / "x_R1_finngen_fatigue_clusters.tsv", sep="\t", index=False, float_format="%.5g")
    m[["cluster_name", "z", "P"]].to_csv(RESULTS / "x_R1_finngen_fatigue_profile.tsv", sep="\t", float_format="%.5g")
    print("supported (H-R):", bool(obs > np.percentile(null, 95)))


if __name__ == "__main__":
    main()
