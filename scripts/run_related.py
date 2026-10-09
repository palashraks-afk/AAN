"""Related conditions (preregistration/PREREG_v1.6_related_conditions.md): UK Biobank fatigue GWAS and long COVID.

Heritability gate first; MAGMA only if z >= 4. Waits for the main overnight run so the CPU is not shared.
"""
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DERIVED = Path("D:/AAN_data/derived")
PANEL = Path("D:/AAN_data/panel")
PY = sys.executable

TRAITS = {
    "mecfs_ukb_donertas": (PANEL / "mecfs_ukb_donertas.h.tsv.gz", 484598),
    "longcovid_strict": (PANEL / "longcovid_strict.h.tsv.gz", 997600),
}


def run(cmd):
    print(" ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], cwd=HERE).returncode == 0


def z_of(name):
    f = DERIVED / "ldsc_results.tsv"
    if not f.exists():
        return None
    t = pd.read_csv(f, sep="\t").drop_duplicates("name", keep="last").set_index("name")
    return float(t.loc[name, "h2_z"]) if name in t.index else None


def main():
    for name, (path, n) in TRAITS.items():
        if z_of(name) is None and path.exists():
            run([PY, "-u", "ldsc.py", "--name", name, "--format", "harmonised", "--gwas", path, "--n", n])
    print("gate:", {n: z_of(n) for n in TRAITS}, flush=True)
    while "overnight run finished" not in (DERIVED / "overnight.log").read_text():
        time.sleep(120)
    for name, (path, n) in TRAITS.items():
        z = z_of(name)
        if z is None or z < 4:
            print(f"{name}: excluded by the heritability gate (z = {z})", flush=True)
            continue
        base = DERIVED / "magma" / name
        if not any((base / f"cluster_B.gsa.out{s}").exists() for s in (".txt", "")):
            run([PY, "-u", "run_magma.py", "--name", name, "--format", "harmonised", "--gwas", path,
                 "--n", n, "--jobs", 6])
        if not (Path("D:/AAN/results") / f"sldsc_{name}.tsv").exists():
            run([PY, "-u", "sldsc_lite.py", "--name", name, "--format", "harmonised", "--gwas", path, "--n", n])


if __name__ == "__main__":
    main()
