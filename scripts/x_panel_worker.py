"""Run MAGMA for panel traits as their downloads complete (file size equals the size in data/panel.tsv).

    python scripts/x_panel_worker.py depression_broad neuroticism als epilepsy ...
"""
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PANEL = Path("D:/AAN_data/panel")
DERIVED = Path("D:/AAN_data/derived")


def done(label):
    return any((DERIVED / "magma" / label / f"cluster_B.gsa.out{s}").exists() for s in (".txt", ""))


def main(labels):
    tab = pd.read_csv(HERE.parent / "data" / "panel.tsv", sep="\t").set_index("label")
    todo = [l for l in labels if not done(l)]
    while todo:
        progressed = False
        for l in list(todo):
            f = PANEL / f"{l}.h.tsv.gz"
            if f.exists() and f.stat().st_size == int(tab.loc[l, "bytes"]):
                subprocess.run([sys.executable, "-u", "run_magma.py", "--name", l, "--format", "harmonised", "--gwas", str(f),
                                "--jobs", "4", "--n", str(int(tab.loc[l, "n_total"]))], cwd=HERE)
                todo.remove(l)
                progressed = True
        if not progressed:
            time.sleep(60)
    print("worker finished", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
