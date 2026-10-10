# HANDOFF: where the project stands and how to continue

Last updated: 2026-10-09 07:35 (Pacific). Read this first, then `TODO_PALASH.md` (what you must do) and `SCORES.md` (honest
grades). Deadline: **Tuesday 2026-10-20, 11:59 p.m. CT.** Aim to submit Monday 19 October.

## 1. Where things stand right now
- Plan, 8 pre-registrations (v1 to v1.6, each tagged and pushed), pipeline code, tests, background research, figure and 3D
  tooling and a report skeleton are all in this repo.
- Heavy computation (MAGMA SNP-to-gene, then cell-type tests):
  | Trait | Role | Status |
  |---|---|---|
  | `decodeme_gwas_1` (ME/CFS primary) | main analysis | MAGMA done, **cell-type table not opened by anyone** |
  | Alzheimer's | positive control (microglia) | done, passes (microglia FDR 3e-7) |
  | Schizophrenia | positive control (excitatory neurons) | done, **not yet checked formally** |
  | Height | negative control (non-neural) | 17 of 22 chromosomes done, last 5 running (slow trait, see section 6) |
  | `decodeme_gwas_2`, female, male, infection subsets | robustness and subtype tests | not done yet |
  | Comparison panel (19 GWAS) | the null for the specificity test | not done yet |
  | FinnGen neglected conditions (4) | extension | heritability gate only; trigeminal neuralgia fails (z = 1.9) |
  | Simplified S-LDSC (own code) | confirmation method | done for the controls, ME/CFS gwas_1 and gwas_2 (files exist, unopened) |
- **Control gate G3: PASS (2026-10-09 19:36).** After that, the first ME/CFS look (T1 only) was taken: see `DISCOVERY.md` section 4. It is
  preliminary: 11 clusters at FDR < 0.05 in model A (excitatory cortical and amygdala neurons, z about 3.6), none in model B.
  **Update 2026-10-10: the full pipeline has finished.** Results and the A to H decisions are in `DISCOVERY.md` section 4. Short version: 6 clusters
  (mainly amygdala excitatory neurons) are specific to ME/CFS against the panel, modest strength; B, D, E not novel or not found; H exploratory.

## 2. Hard rules (do not break these)
1. Do not open or interpret ME/CFS cluster-level results until `check_controls.py` prints `G3 overall: PASS`. If it fails,
   debug the pipeline, not the ME/CFS result. (`.gitignore` keeps the ME/CFS cell-type tables out of git until then.)
2. Every change to a pre-registered threshold or definition goes in `preregistration/DEVIATIONS.md` with date and reason, and
   only if it was decided before the relevant result was seen.
3. Never write "first", "never done", "groundbreaking", "cure" or "treatment" in the application. Write "to my knowledge not
   found in my search of [date]". Findings are hypotheses for researchers, not advice.
4. The prize needs your own written work. `paper/DRAFT_report.md` is scaffolding: rewrite everything in your own words.
5. Don't commit data or large outputs. Everything big lives in `AAN_data`.
6. **Before any `git add -A`, run `git status` and `git diff --stat`.** A commit made from a stale working tree once deleted 14
   files (see `docs/10_self_review.md`, item 10). Several scripts and notes were written by different sessions.

## 3. The computer must stay awake
**This has cost about 32 hours so far.** Jobs were killed by sleep at 01:26 and 23:20 on 8 Oct and at 07:54 on 9 Oct (the app's own
keep-awake only lasts about 5 minutes after a session goes idle). Two protections now exist:
1. `scripts/keep_awake.ps1`, a detached process that holds Windows awake until `pipeline.log` says "pipeline finished". Start it with
   `powershell -NoProfile -WindowStyle Hidden -File D:\AAN\scripts\keep_awake.ps1`. It changes no settings; end it by stopping the process.
2. **You should also set sleep to "never" while plugged in** (Windows: Settings, System, Power) and keep the lid open: a closed lid sleeps
   whatever any program asks.

Jobs have been killed twice by the computer sleeping (01:26 on 8 Oct, 23:20 on 8 Oct). Keep it plugged in, lid open, sleep set
to "never" while plugged in (Windows: Settings, System, Power). Everything resumes: `run_magma.py` skips chromosomes that are
already finished, and the other scripts skip finished traits.

## 3b. Two machines are running this project (Windows and Mac)
Git history shows a second session rebuilt and ran the pipeline on a Mac (commits of 2026-10-08 and 2026-10-09). Data and results live
outside git, so the two machines do **not** share outputs. State as last recorded:
| | Windows (this repo's `D:\AAN_data`) | Mac (per `docs/10_self_review.md`) |
|---|---|---|
| ME/CFS gwas_1 MAGMA | done | done |
| other five ME/CFS GWAS | not done | done |
| Alzheimer's | done | done |
| Schizophrenia | done | still running when last noted |
| Height | 17 of 22 chromosomes, resuming | still running when last noted |
Rules so the work does not clash:
1. **Pick one machine as authoritative for the final results.** Whichever passes the control gate G3 first, with the fewest steps left.
   Do not mix result files from the two machines in one analysis.
2. Both run the same code (`run_everything.py`), the same p-value floor (1e-300) and the same pre-registrations, so the two should agree. If both
   finish, compare `genes.genes.out` for `decodeme_gwas_1` (expect near-identical z-scores); a mismatch is itself a finding to log.
3. Before any push: `git pull --rebase origin main`, then `git status` and `git diff --stat`. Pushes were rejected once because the
   other machine had pushed first.
4. Tag the machine in `results/RUN_REPORT.md` (the report shows the date; add the machine name by hand).

## 4. How to run everything that is left: one command
```
cd D:\AAN\scripts
python run_everything.py --dry-run      # shows every step and whether it is already DONE or still TODO
python run_everything.py                # does everything that is left, in order, then writes results/RUN_REPORT.md
```
It is safe to start again at any time (after sleep, a crash, or a new session): a step whose outputs exist is skipped and MAGMA
resumes by chromosome. It never commits to git. Progress: `D:\AAN_data\derived\pipeline.log`.

Order inside the runner:
1. ME/CFS primary and the three controls (schizophrenia, Alzheimer's, height).
2. **Control gate G3** (`check_controls.py`). If it prints FAIL the run stops and writes `results/G3_FAILED.txt`; nothing that reads
   ME/CFS cell-type results runs. If PASS it continues.
3. All remaining compute: the other five ME/CFS GWAS, heritability for all, the 19-trait comparison panel, the FinnGen neglected
   conditions, long COVID and the UK Biobank fatigue GWAS, simplified S-LDSC, the method benchmark (v1.4).
4. Discovery-layer MAGMA runs (groups, whole-body cell types, regions, Reactome).
5. Analysis of ME/CFS (only after G3 PASS): calibration for all six GWAS, conditioning (model C), locus and chromosome drop,
   permutation null, per-locus contributions, independent signals, discovery layers, gene ranking and benchmark on migraine,
   rheumatoid arthritis and IBD, effect direction, 3D protein models of the top genes, 3D cell landscape, figures (including the 3D
   brain), and the evidence ledger.
6. `results/RUN_REPORT.md` lists every step as done, skipped or failed.

To run just some steps: `python run_everything.py --only calibrate_gwas_1 figures` (names as in `--dry-run`).

After it finishes, the manual part: read `results/EVIDENCE_LEDGER.md`, fill `DISCOVERY.md` section 4 and the report Results from the
real numbers, then write the abstract (under 300 words) last, in your own words.

## 5. Where everything is
| Thing | Path | Notes |
|---|---|---|
| Repo | `D:\AAN` | this repository (github.com/palashraks-afk/AAN) |
| Data and big outputs | `D:\AAN_data` | about 60 GB, never committed |
| Python environment | `D:\AAN_data\venv` | rebuild from `requirements.txt` |
| Reference genotypes by chromosome | `C:\AAN_ref\by_chr` | on the SSD, made by `split_reference.py` |
| MAGMA | `D:\AAN_data\tools\magma\magma.exe` | v1.10 Windows build |
| MAGMA results per trait | `D:\AAN_data\derived\magma\<trait>` | `cluster_A/B.gsa.out.txt` are the cell-type tests |
| Run log | `D:\AAN_data\derived\overnight.log` and `overnight_stdout.log` | check progress here |
| Paths are hard-coded as `D:/AAN_data`, `C:/AAN_ref`, `D:/AAN` | | on another machine, see `docs/09_new_machine_setup.md` |

### Script index
| Script | Does |
|---|---|
| `download_decodeme.py`, `download_panel.py`, `make_panel.py` | fetch and list the data (MD5 or integrity checked) |
| `build_rsid_map.py`, `build_specificity.py`, `build_region_covars.py`, `build_discovery_covars.py`, `split_reference.py` | inputs for the analysis |
| `ldsc.py` | own LD score regression (heritability gate); tests in `tests/` |
| **`run_everything.py`** | **the one command that runs every remaining step in order, with the control gate** (older `run_overnight.py` does only the compute part) |
| `run_magma.py`, `gene_property.py`, `run_panel.py`, `run_extension.py`, `run_related.py`, `run_layers.py`, `method_benchmark.py` | the MAGMA runs, discovery layers and the v1.4 method benchmark |
| `sldsc_lite.py`, `run_sldsc_batch.py` | simplified partitioned heritability (confirmation method C1) |
| `check_controls.py`, `check_split_reference.py` | the control gate G3; split-reference equivalence check |
| `calibrate.py`, `analyze_layers.py`, `independent_signals.py`, `confirm_robustness.py`, `make_evidence_ledger.py` | results, discovery layers, confirmation, evidence ledger |
| `rank_targets.py`, `get_known_targets.py`, `smr_lite.py` | gene ranking, benchmark targets, effect direction |
| `make_figures.py`, `brain3d.py`, `landscape3d.py`, `structure3d.py`, `background_stats.py` | figures, 3D brain, 3D cell landscape, 3D protein models, background numbers |

## 6. Known problems and what fixed them
- DecodeME files are raw REGENIE output; only variants in `gwas_qced.var.gz` are used (logged in `DEVIATIONS.md`).
- The standard `ldsc` package will not build on Windows (needs pysam): own implementation in `ldsc.py`.
- MAGMA with the full reference was slow; the per-chromosome reference is about 3 times faster and gives identical results.
- Windows MAGMA writes `*.gsa.out.txt`; the scripts look for both names.
- **Height is very slow in MAGMA** (about 8 hours): 26,000 SNPs have p below 1e-30 and some underflow to 5e-324. Capping p-values
  only halves the time, so the pre-registered setup was kept.
- Background jobs die if the computer sleeps or the session ends; see section 3.
- One failed trait no longer stops the panel run. Both harmonised file layouts are read.

## 7. Prior-art re-search still to do
Run on **17 October**: PubMed, bioRxiv and medRxiv for "ME/CFS" with cell type, single-nucleus, LDSC, specificity, brainstem,
hypothalamus, plus whether the DecodeME paper is now published with new analyses, plus the larger ME/CFS meta-analysis
(46,450 cases, newsletter report only). Record the date and terms in `docs/02_prior_art_and_sources.md`.
