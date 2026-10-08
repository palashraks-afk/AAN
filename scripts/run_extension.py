"""Neglected-condition extension (preregistration/PREREG_v1.1_extension.md).

For each FinnGen condition: heritability gate first (cheap), MAGMA only if it passes (z >= 4).
MAGMA waits until the main overnight run has finished so the two don't fight over the CPU.
"""
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DERIVED = Path("D:/AAN_data/derived")
DATA = Path("D:/AAN_data/panel/neglected")
PY = sys.executable
N_TOTAL = 480_000

CONDITIONS = {
    "trigeminal_neuralgia": "G6_TRINEU",
    "cluster_headache": "G6_CLUSTHEADACHE_WIDE",
    "menieres": "H8_MENIERE",
    "dystonia": "G6_DYSTON",
}


def run(cmd):
    print(" ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], cwd=HERE).returncode == 0


def ldsc_z(name):
    f = DERIVED / "ldsc_results.tsv"
    if not f.exists():
        return None
    t = pd.read_csv(f, sep="\t").drop_duplicates("name", keep="last").set_index("name")
    return float(t.loc[name, "h2_z"]) if name in t.index else None


def main():
    for name, endpoint in CONDITIONS.items():
        path = DATA / f"finngen_R13_{endpoint}.gz"
        if ldsc_z(name) is None:
            run([PY, "-u", "ldsc.py", "--name", name, "--format", "finngen", "--gwas", path, "--n", N_TOTAL])
    print("heritability gate:", {n: ldsc_z(n) for n in CONDITIONS}, flush=True)

    while "overnight run finished" not in (DERIVED / "overnight.log").read_text():
        time.sleep(120)

    for name, endpoint in CONDITIONS.items():
        z = ldsc_z(name)
        if z is None or z < 4:
            print(f"{name}: excluded by the heritability gate (z = {z})", flush=True)
            continue
        base = DERIVED / "magma" / name
        if not any((base / f"cluster_B.gsa.out{s}").exists() for s in (".txt", "")):
            run([PY, "-u", "run_magma.py", "--name", name, "--format", "finngen",
                 "--gwas", DATA / f"finngen_R13_{endpoint}.gz", "--n", N_TOTAL, "--jobs", 6])


if __name__ == "__main__":
    main()
