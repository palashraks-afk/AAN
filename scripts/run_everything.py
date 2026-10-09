"""Run every step of the project that is not finished yet, in order, with the control gate built in.

    python scripts/run_everything.py              # do everything that is left
    python scripts/run_everything.py --dry-run    # print the plan and what is already done
    python scripts/run_everything.py --only calibrate independent_signals

Safe to start again at any time: a step whose outputs exist is skipped, and MAGMA resumes by chromosome.
Order: compute for the ME/CFS GWAS and the three controls -> control gate G3 -> all remaining compute -> discovery layers ->
analysis of ME/CFS. The analysis steps read ME/CFS cell-type results, so they only run after G3 prints PASS (rule 1 in
HANDOFF.md). If G3 fails the run stops and writes results/G3_FAILED.txt. Nothing here commits to git.
Progress goes to D:/AAN_data/derived/pipeline.log and the final summary to results/RUN_REPORT.md.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DERIVED = Path("D:/AAN_data/derived")
DECODEME = Path("D:/AAN_data/decodeme")
PANEL_DIR = Path("D:/AAN_data/panel")
RES = ROOT / "results"
FIG = ROOT / "figures"
PY = sys.executable
LOG = DERIVED / "pipeline.log"
OLD_LOG = DERIVED / "overnight.log"

CONTROLS = ["schizophrenia", "alzheimer", "height"]
DECODEME_GWAS = ["gwas_1", "gwas_2", "gwas_1_female", "gwas_1_infectious_onset", "gwas_1_non_infectious_onset", "gwas_1_male"]
BENCHMARK_DISEASES = ["migraine", "rheumatoid_arthritis", "ibd"]
INDEPENDENT_SIGNAL_TRAITS = ["decodeme_gwas_1", "schizophrenia", "alzheimer", "rheumatoid_arthritis"]
TARGET = "decodeme_gwas_1"

results = {}          # step name -> "done", "skipped (reason)", "failed"


def say(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def magma_done(name):
    base = DERIVED / "magma" / name
    return any((base / f"cluster_B.gsa.out{s}").exists() for s in (".txt", ""))


def ldsc_done(name):
    f = DERIVED / "ldsc_results.tsv"
    if not f.exists():
        return False
    return name in pd.read_csv(f, sep="\t")["name"].tolist()


def exists(*paths):
    return all(Path(p).exists() for p in paths)


def downloaded_panel_labels():
    panel = pd.read_csv(ROOT / "data" / "panel.tsv", sep="\t")
    labels = []
    for label in panel["label"]:
        if label == "schizophrenia":
            continue
        if (PANEL_DIR / f"{label}.h.tsv.gz").exists():
            labels.append(label)
    return labels


def panel_done():
    return all(magma_done(l) and ldsc_done(l) for l in downloaded_panel_labels())


def run(cmd, cwd=HERE):
    say("RUN " + " ".join(str(c) for c in cmd))
    code = subprocess.run([str(c) for c in cmd], cwd=cwd).returncode
    if code != 0:
        say(f"  step exited with code {code}")
    return code == 0


def wait_for_magma(name):
    """If an earlier run left MAGMA jobs for this trait running, let them finish instead of starting duplicates."""
    marker = f"magma\\{name}\\".lower()
    while True:
        ps = subprocess.run(["powershell", "-NoProfile", "-Command",
                             "Get-CimInstance Win32_Process -Filter \"Name='magma.exe'\" | "
                             "Select-Object -ExpandProperty CommandLine"], capture_output=True, text=True).stdout.lower()
        if marker not in ps.replace("/", "\\"):
            return
        say(f"waiting for running MAGMA jobs of {name}")
        time.sleep(60)


def candidate_symbols(n=6):
    f = RES / f"candidate_genes_{TARGET}.tsv"
    if not f.exists():
        return []
    return pd.read_csv(f, sep="\t")["symbol"].dropna().head(n).tolist()


def build_steps():
    S = []

    def add(name, cmd, done, gated=False, stage="compute"):
        S.append({"name": name, "cmd": cmd, "done": done, "gated": gated, "stage": stage})

    # ---- stage 1: ME/CFS primary and the controls
    add("magma_decodeme_gwas_1",
        [PY, "-u", "run_magma.py", "--name", "decodeme_gwas_1", "--format", "decodeme",
         "--gwas", DECODEME / "gwas_1.regenie.gz", "--jobs", 6, "--reuse"],
        lambda: magma_done("decodeme_gwas_1"))
    add("controls_magma", [PY, "-u", "run_panel.py", "--only"] + CONTROLS + ["--skip-ldsc"],
        lambda: all(magma_done(c) for c in CONTROLS))
    # the gate is a special step handled in main()

    # ---- stage 2: everything else that needs the CPU
    for stem in DECODEME_GWAS[1:]:
        add(f"magma_decodeme_{stem}",
            [PY, "-u", "run_magma.py", "--name", f"decodeme_{stem}", "--format", "decodeme",
             "--gwas", DECODEME / f"{stem}.regenie.gz", "--jobs", 6, "--reuse"],
            (lambda s=stem: magma_done(f"decodeme_{s}")))
    for stem in DECODEME_GWAS:
        add(f"ldsc_decodeme_{stem}",
            [PY, "-u", "ldsc.py", "--name", f"decodeme_{stem}", "--format", "decodeme",
             "--gwas", DECODEME / f"{stem}.regenie.gz"],
            (lambda s=stem: ldsc_done(f"decodeme_{s}")))
    add("panel_all", [PY, "-u", "run_panel.py"], panel_done)
    add("extension_neglected_conditions", [PY, "-u", "run_extension.py"],
        lambda: all(ldsc_done(n) for n in ("trigeminal_neuralgia", "cluster_headache", "menieres", "dystonia")))
    add("related_conditions", [PY, "-u", "run_related.py"], lambda: ldsc_done("longcovid_strict") and ldsc_done("mecfs_ukb_donertas"))
    add("sldsc_batch", [PY, "-u", "run_sldsc_batch.py"], lambda: exists(RES / f"sldsc_{TARGET}.tsv", RES / "sldsc_height.tsv"))
    add("method_benchmark", [PY, "-u", "method_benchmark.py"], lambda: exists(RES / "method_benchmark.md"), stage="analysis-nogate")

    # ---- stage 3: discovery-layer MAGMA runs (no ME/CFS result is read)
    add("region_inputs", [PY, "-u", "build_region_covars.py"], lambda: exists(DERIVED / "gene_covar_regions.txt"))
    add("layers_magma", [PY, "-u", "run_layers.py"],
        lambda: exists(DERIVED / "magma" / TARGET / "reactome_C.gsa.out.txt", DERIVED / "magma" / TARGET / "regions_A.gsa.out.txt"))

    # ---- stage 4: analysis of ME/CFS (after G3)
    for stem in DECODEME_GWAS:
        add(f"calibrate_{stem}", [PY, "-u", "calibrate.py", "--target", f"decodeme_{stem}"],
            (lambda s=stem: exists(RES / f"clusters_decodeme_{s}.tsv")), gated=True, stage="analysis")
    add("conditioning_inputs", [PY, "-u", "make_conditioning_covar.py"], lambda: exists(DERIVED / "gene_covar_C.txt"),
        gated=True, stage="analysis")
    add("model_C", [PY, "-u", "gene_property.py", "--name", TARGET, "--model", "C"],
        lambda: any((DERIVED / "magma" / TARGET / f"cluster_C.gsa.out{s}").exists() for s in (".txt", "")),
        gated=True, stage="analysis")
    add("locus_drop", [PY, "-u", "confirm_robustness.py", "--name", TARGET, "--step", "locus"],
        lambda: any((DERIVED / "magma" / TARGET / f"cluster_loco_chr22.gsa.out{s}").exists() for s in (".txt", "")),
        gated=True, stage="analysis")
    add("permutation_null", [PY, "-u", "confirm_robustness.py", "--name", TARGET, "--step", "permute"],
        lambda: exists(RES / f"permutation_{TARGET}.tsv"), gated=True, stage="analysis")
    add("per_locus", [PY, "-u", "confirm_robustness.py", "--name", TARGET, "--step", "perlocus"],
        lambda: exists(RES / f"perlocus_{TARGET}.tsv"), gated=True, stage="analysis")
    for t in INDEPENDENT_SIGNAL_TRAITS:
        add(f"independent_signals_{t}", [PY, "-u", "independent_signals.py", "--name", t],
            (lambda n=t: exists(RES / f"independent_signals_{n}.tsv")), gated=(t == TARGET), stage="analysis")
    add("layers_analysis", [PY, "-u", "analyze_layers.py"], lambda: exists(RES / f"groups_A_{TARGET}.tsv", RES / "comorbidity_map.tsv"),
        gated=True, stage="analysis")
    add("rank_candidates", [PY, "-u", "rank_targets.py", "--trait", TARGET],
        lambda: exists(RES / f"candidate_genes_{TARGET}.tsv"), gated=True, stage="analysis")
    for d in BENCHMARK_DISEASES:
        add(f"benchmark_{d}", [PY, "-u", "rank_targets.py", "--trait", d, "--disease", d],
            (lambda n=d: (RES / "benchmark.tsv").exists() and n in (RES / "benchmark.tsv").read_text()), stage="analysis-nogate")
    add("effect_direction", [PY, "-u", "smr_lite.py", "--candidates", RES / f"candidate_genes_{TARGET}.tsv", "--top", 200],
        lambda: exists(RES / f"smr_lite_{TARGET}.tsv"), gated=True, stage="analysis")
    add("structures_3d", lambda: [PY, "-u", "structure3d.py"] + candidate_symbols(),
        lambda: bool(candidate_symbols()) and all(exists(FIG / "structures" / f"{s}.png") for s in candidate_symbols()),
        gated=True, stage="analysis")
    add("landscape_3d", [PY, "-u", "landscape3d.py", "--trait", TARGET], lambda: exists(FIG / "fig8_landscape_decodeme_gwas_1.png"),
        gated=True, stage="analysis")
    add("figures", [PY, "-u", "make_figures.py"], lambda: exists(FIG / "fig2_enrichment.png", FIG / "fig7_brain3d.png"),
        gated=True, stage="analysis")
    add("evidence_ledger", [PY, "-u", "make_evidence_ledger.py"], lambda: exists(RES / "EVIDENCE_LEDGER.md"),
        gated=True, stage="analysis")
    return S


def gate_g3():
    """Run check_controls.py and return True only if it prints 'G3 overall: PASS'."""
    say("RUN control gate G3 (check_controls.py)")
    out = subprocess.run([PY, "check_controls.py"], cwd=HERE, capture_output=True, text=True)
    text = out.stdout + out.stderr
    (RES / "G3_check.txt").write_text(text, encoding="utf-8")
    passed = "G3 overall: PASS" in text
    say("G3 " + ("PASS" if passed else "FAIL"))
    if not passed:
        (RES / "G3_FAILED.txt").write_text("The control gate failed. Do not read the ME/CFS results; debug the pipeline.\n\n" + text,
                                           encoding="utf-8")
    return passed


def write_report(gate):
    lines = ["# Run report", "", f"Generated {time.strftime('%Y-%m-%d %H:%M:%S')} by scripts/run_everything.py.", "",
             f"Control gate G3: **{gate}**", "", "| step | status |", "|---|---|"]
    for name, status in results.items():
        lines.append(f"| {name} | {status} |")
    failed = [n for n, s in results.items() if s.startswith("failed")]
    lines += ["", f"{len(failed)} step(s) failed. Start the script again to retry them; finished steps are skipped." if failed
              else "No step failed.", ""]
    lines.append("Next: read `results/EVIDENCE_LEDGER.md`, then fill `DISCOVERY.md` section 4 and the report Results from the numbers.")
    (RES / "RUN_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", nargs="*", help="run just these steps (names as in --dry-run)")
    args = ap.parse_args()
    RES.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    steps = build_steps()

    if args.dry_run:
        for s in steps:
            print(f"{'DONE ' if s['done']() else 'TODO '} {s['stage']:<16} {'[after G3] ' if s['gated'] else '           '}{s['name']}")
        return 0

    say("pipeline started")
    gate = {"value": None}            # None = not checked yet, True = passed, False = failed

    def ensure_gate():
        if gate["value"] is None and all(magma_done(c) for c in CONTROLS):
            gate["value"] = gate_g3()
        return gate["value"]

    def release_old_log_phrase():
        # run_extension.py and run_related.py wait for this phrase in the old log before they start MAGMA
        if "overnight run finished" not in OLD_LOG.read_text():
            with open(OLD_LOG, "a") as fh:
                fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + "overnight run finished" + chr(10))

    for s in steps:
        if args.only and s["name"] not in args.only:
            continue
        if s["gated"]:
            g = ensure_gate()
            if g is None:
                results[s["name"]] = "skipped (controls not finished)"
                continue
            if g is False:
                results[s["name"]] = "skipped (gate failed)"
                continue
        if s["done"]():
            results[s["name"]] = "done (already)"
        else:
            if s["name"] == "controls_magma":
                for c in CONTROLS:
                    wait_for_magma(c)
            if s["name"].startswith(("extension_", "related_")):
                release_old_log_phrase()
            cmd = s["cmd"]() if callable(s["cmd"]) else s["cmd"]
            ok = run(cmd)
            results[s["name"]] = ("done" if s["done"]() else "done (no output check)") if ok else "failed"
        if s["name"] == "controls_magma" and ensure_gate() is False:
            say("stopping: the control gate failed; fix the pipeline, then start again")
            write_report("FAIL")
            return 1
    write_report({None: "not reached (controls not finished)", True: "PASS", False: "FAIL"}[gate["value"]])
    say("pipeline finished")
    return 0


if __name__ == "__main__":
    sys.exit(main())
