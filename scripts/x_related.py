"""NEXT_20 M4 (PREREG v1.15): is the six-cluster amygdala/cortical pattern present in related conditions?

    python scripts/x_related.py
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
REL = Path("D:/AAN_data/panel/related")
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]
N_TOTAL = 480_000
JOBS = {  # name: (format, file)
    "finngen_fibromyalgia": ("finngen", REL / "finngen_R13_M13_FIBROMYALGIA.gz"),
    "finngen_pain": ("finngen", REL / "finngen_R13_PAIN.gz"),
    "ibs": ("existing", None),
}


def run(cmd):
    print(" ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], cwd=HERE).returncode == 0


def gate(name, fmt, path):
    f = DERIVED / "ldsc_results.tsv"
    have = f.exists() and name in pd.read_csv(f, sep="\t")["name"].values
    if not have:
        run([sys.executable, "-u", "ldsc.py", "--name", name, "--format", fmt, "--gwas", path, "--n", N_TOTAL])
    t = pd.read_csv(f, sep="\t").drop_duplicates("name", keep="last").set_index("name")
    return float(t.loc[name, "h2_z"])


def main():
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    rows = []
    for name, (fmt, path) in JOBS.items():
        if fmt == "existing":
            z_h2 = None
        else:
            z_h2 = gate(name, fmt, path)
            print(name, "heritability z", z_h2)
            if z_h2 < 4:
                rows.append({"condition": name, "h2_z": z_h2, "note": "excluded by the heritability gate"})
                continue
            if not any((DERIVED / "magma" / name / f"cluster_A.gsa.out{s}").exists() for s in (".txt", "")):
                run([sys.executable, "-u", "run_magma.py", "--name", name, "--format", fmt, "--gwas", path, "--n", N_TOTAL, "--jobs", 8])
        m = calibrate.load_model(name, "A")
        m["z"] = calibrate.z_from_p(m["P"].to_numpy())
        m = m.join(ann[["cluster_name"]])
        six = m[m["cluster_name"].isin(SIX)]
        obs = six["z"].mean()
        rng = np.random.default_rng(20261010)
        null = np.array([rng.choice(m["z"].to_numpy(), 6, replace=False).mean() for _ in range(10000)])
        p = (1 + (null >= obs).sum()) / 10001
        rows.append({"condition": name, "h2_z": z_h2, "mean_z_six": round(obs, 3), "null_95th": round(np.percentile(null, 95), 3), "p": round(p, 4),
                     "n_positive_of_6": int((six["z"] > 0).sum()), "supported_p_lt_0.0167": bool(p < 0.05 / 3),
                     **{f"z_{r.cluster_name}": round(r.z, 2) for r in six.itertuples()}})
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "x_M4_related_conditions.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
