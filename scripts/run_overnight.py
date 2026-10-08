"""Run the whole heavy part of the pipeline in a sensible order. Safe to restart: finished steps are skipped.

Order: ME/CFS gwas_1 first, then the three control traits (so G3 can be checked early), then the other
five ME/CFS GWAS, then the heritability gate for all of them, then the rest of the comparison panel.

    python scripts/run_overnight.py
"""
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DERIVED = Path("D:/AAN_data/derived")
DECODEME = Path("D:/AAN_data/decodeme")
PY = sys.executable
LOG = DERIVED / "overnight.log"

DECODEME_GWAS = ["gwas_1", "gwas_2", "gwas_1_female", "gwas_1_infectious_onset",
                 "gwas_1_non_infectious_onset", "gwas_1_male"]
CONTROLS = ["schizophrenia", "alzheimer", "height"]


def say(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def magma_done(name):
    base = DERIVED / "magma" / name
    return any((base / f"cluster_B.gsa.out{s}").exists() for s in (".txt", ""))


def ldsc_done(name):
    f = DERIVED / "ldsc_results.tsv"
    return f.exists() and f"\t{name}" in f.read_text() or (f.exists() and f.read_text().rstrip().endswith(name))


def run(cmd):
    say("RUN " + " ".join(str(c) for c in cmd))
    result = subprocess.run([str(c) for c in cmd], cwd=HERE)
    if result.returncode != 0:
        say(f"FAILED (exit {result.returncode}), continuing with the next step")
    return result.returncode == 0


def decodeme_magma(stem):
    name = f"decodeme_{stem}"
    if not magma_done(name):
        run([PY, "-u", "run_magma.py", "--name", name, "--format", "decodeme",
             "--gwas", DECODEME / f"{stem}.regenie.gz", "--jobs", 6])


def decodeme_ldsc(stem):
    name = f"decodeme_{stem}"
    if not ldsc_done(name):
        run([PY, "-u", "ldsc.py", "--name", name, "--format", "decodeme",
             "--gwas", DECODEME / f"{stem}.regenie.gz"])


def main():
    say("overnight run started")
    decodeme_magma("gwas_1")
    run([PY, "-u", "run_panel.py", "--only"] + CONTROLS + ["--skip-ldsc"])
    for stem in DECODEME_GWAS[1:]:
        decodeme_magma(stem)
    for stem in DECODEME_GWAS:
        decodeme_ldsc(stem)
    run([PY, "-u", "run_panel.py"])
    say("overnight run finished")


if __name__ == "__main__":
    main()
