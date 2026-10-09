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
- **Nothing about ME/CFS cell types has been read.** The rule: the controls must pass first (`check_controls.py` must print
  `G3 overall: PASS`). That check takes seconds once height finishes.

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
Jobs have been killed twice by the computer sleeping (01:26 on 8 Oct, 23:20 on 8 Oct). Keep it plugged in, lid open, sleep set
to "never" while plugged in (Windows: Settings, System, Power). Everything resumes: `run_magma.py` skips chromosomes that are
already finished, and the other scripts skip finished traits.

## 4. What to do next, in order
Start (or restart after any interruption) with:
```
cd D:\AAN\scripts
python run_overnight.py            # finishes height, then the other ME/CFS GWAS, LDSC for all, then the panel
python run_extension.py            # FinnGen neglected conditions, waits for the above
python run_related.py              # long COVID and UK Biobank fatigue, outside the null panel
```
Then, in this order:
1. `python check_controls.py` needs `G3 overall: PASS` (schizophrenia, Alzheimer's, height). If FAIL: stop and debug.
2. `python calibrate.py --target decodeme_gwas_1`, the ME/CFS cell-type table with specificity calibration (T1 and T2).
3. `python make_conditioning_covar.py` then `python gene_property.py --name decodeme_gwas_1 --model C` (T3; needs the
   depression, BMI and insomnia panel traits finished).
4. Confirmation: `python confirm_robustness.py --name decodeme_gwas_1 --step locus`, then `--step permute`; S-LDSC results
   already exist; then `python make_evidence_ledger.py` (writes `results/EVIDENCE_LEDGER.md`).
5. Discovery layers: `python build_region_covars.py`, `python run_layers.py`, `python analyze_layers.py`,
   `python independent_signals.py --name decodeme_gwas_1`, `python landscape3d.py --trait decodeme_gwas_1`.
6. Translational: `python rank_targets.py --trait decodeme_gwas_1`; benchmark with `--trait migraine --disease migraine`,
   `rheumatoid_arthritis`, `ibd`; then `python smr_lite.py --candidates results/candidate_genes_decodeme_gwas_1.tsv`;
   `python structure3d.py GENE1 GENE2 ...` for 3D protein models.
7. Figures: `python make_figures.py` (3D brain: `figures/fig7_brain3d.png` and `.html`).
8. Fill `DISCOVERY.md` section 4 and `paper/DRAFT_report.md` Results from the real numbers, then write the abstract (under 300
   words) last, in your own words.

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
| `run_magma.py`, `gene_property.py`, `run_panel.py`, `run_overnight.py`, `run_extension.py`, `run_related.py`, `run_layers.py` | the MAGMA runs and orchestration |
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
