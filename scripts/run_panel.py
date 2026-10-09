"""Run LDSC and MAGMA for every downloaded panel trait, one after another.

    python scripts/run_panel.py                 # everything that is downloaded and not yet done
    python scripts/run_panel.py --only height alzheimer schizophrenia

Skips a trait if its MAGMA output already exists, so it is safe to restart after an interruption.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PANEL = Path("D:/AAN_data/panel")
DERIVED = Path("D:/AAN_data/derived")
PY = sys.executable


def run(cmd):
    """Run one step; a failure is logged and the loop carries on with the next trait."""
    print(" ".join(str(c) for c in cmd), flush=True)
    result = subprocess.run([str(c) for c in cmd], cwd=HERE)
    if result.returncode != 0:
        print(f"STEP FAILED (exit {result.returncode}): {cmd[1:4]}", flush=True)
    return result.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--skip-ldsc", action="store_true")
    args = ap.parse_args()

    panel = pd.read_csv(HERE.parent / "data" / "panel.tsv", sep="\t")
    jobs = []
    for _, row in panel.iterrows():
        if row["label"] == "schizophrenia":       # replaced by the PGC3 release below
            continue
        jobs.append((row["label"], PANEL / f"{row['label']}.h.tsv.gz", "harmonised", int(row["n_total"])))
    jobs.append(("schizophrenia", PANEL / "PGC3_SCZ_eur.tsv.gz", "pgc3", None))

    for label, path, fmt, n in jobs:
        if args.only and label not in args.only:
            continue
        if not path.exists():
            print("not downloaded yet:", label)
            continue
        extra = ["--n", n] if n else []
        if not args.skip_ldsc:
            done = (DERIVED / "ldsc_results.tsv").exists() and label in (DERIVED / "ldsc_results.tsv").read_text()
            if not done:
                run([PY, "ldsc.py", "--name", label, "--format", fmt, "--gwas", path] + extra)
        if not any((DERIVED / "magma" / label / f"cluster_B.gsa.out{s}").exists() for s in (".txt", "")):
            run([PY, "run_magma.py", "--name", label, "--format", fmt, "--gwas", path, "--jobs", 6] + extra)


if __name__ == "__main__":
    main()
