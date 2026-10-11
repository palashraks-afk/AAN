"""PREREG v1.14 W10: rerun the ME/CFS gene analysis and model A with other SNP-to-gene windows (35/10 kb is the main run)."""
import shutil
import subprocess
import sys
from pathlib import Path

DERIVED = Path("D:/AAN_data/derived")
MAGMA = "D:/AAN_data/tools/magma/magma.exe"
GENELOC = "D:/AAN_data/tools/ref/NCBI38.gene.loc"
PY = sys.executable
SRC = DERIVED / "magma" / "decodeme_gwas_1"

for tag, win in (("w0", "0,0"), ("w50", "50,50")):
    name = f"decodeme_gwas_1_{tag}"
    wd = DERIVED / "magma" / name
    wd.mkdir(parents=True, exist_ok=True)
    for f in ("pval.txt", "snploc.txt"):
        shutil.copy(SRC / f, wd / f)
    subprocess.run([MAGMA, "--annotate", f"window={win}", "--snp-loc", str(wd / "snploc.txt"),
                    "--gene-loc", GENELOC, "--out", str(wd / "annot")], check=True)
    subprocess.run([PY, "-u", "run_magma.py", "--name", name, "--gwas", "unused", "--format", "decodeme", "--reuse"],
                   check=True, cwd="D:/AAN/scripts")
print("W10 done")
