"""Run the simplified S-LDSC (C1 and the M1 benchmark) for the controls and the ME/CFS GWAS, one after another.

Skips traits whose result file already exists. Intended to be started detached so it survives a session restart.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PANEL = Path("D:/AAN_data/panel")
DECODEME = Path("D:/AAN_data/decodeme")
RES = Path("D:/AAN/results")
PY = sys.executable

JOBS = [
    ("alzheimer", "harmonised", PANEL / "alzheimer.h.tsv.gz", 487511),
    ("height", "harmonised", PANEL / "height.h.tsv.gz", 2200007),
    ("decodeme_gwas_1", "decodeme", DECODEME / "gwas_1.regenie.gz", None),
    ("decodeme_gwas_2", "decodeme", DECODEME / "gwas_2.regenie.gz", None),
    ("rheumatoid_arthritis", "harmonised", PANEL / "rheumatoid_arthritis.h.tsv.gz", 79799),
    ("ibd", "harmonised", PANEL / "ibd.h.tsv.gz", 59957),
]

for name, fmt, path, n in JOBS:
    if (RES / f"sldsc_{name}.tsv").exists():
        continue
    cmd = [PY, "-u", "sldsc_lite.py", "--name", name, "--format", fmt, "--gwas", str(path)]
    if n:
        cmd += ["--n", str(n)]
    print("RUN", name, flush=True)
    subprocess.run(cmd, cwd=HERE)
