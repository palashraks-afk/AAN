# HANDOFF: where the project stands and how to continue on another device

Written 2026-10-08 03:20 (Pacific). Read this first. Deadline: **Tuesday 2026-10-20, 11:59 p.m. CT**. Aim to submit Mon 19 Oct.

## 1. One-paragraph state
Plan, pre-registrations (v1 to v1.4, tagged and pushed), pipeline code, background research, figures tooling and a report
skeleton are in this repo. The heavy computation (MAGMA SNP-to-gene for each GWAS, then cell-type tests) is partly done:
ME/CFS primary GWAS (`decodeme_gwas_1`) and Alzheimer's (control) are finished. **No ME/CFS cell-type result has been opened
by anyone**; the rule is that the three control traits must pass first (`scripts/check_controls.py`), and height and
schizophrenia have not finished their MAGMA run yet. After that gate the order of work is in section 4.

## 2. Hard rules (do not break these)
1. Do not open or interpret ME/CFS cluster-level results until `check_controls.py` prints `G3 overall: PASS`. If it fails,
   debug the pipeline, not the ME/CFS result.
2. Every change to a pre-registered threshold or definition is logged in `preregistration/DEVIATIONS.md` with date and reason,
   and only if made before the relevant result was seen.
3. Never write "first", "never done", "cure", "treatment" in the application. Say "to my knowledge, not found in my search
   of [date]". Findings are hypotheses for researchers.
4. The prize needs your own written work. The drafts in `paper/` are scaffolding: rewrite everything in your own words.
5. Don't commit data or large outputs. Everything big lives in `AAN_data` (see section 5).

## 3. What is done and where
| Item | Status | Where |
|---|---|---|
| Plan, grades, prize requirements, schedule | done | `PLAN.md`, `docs/01` to `docs/08` |
| Pre-registrations v1, v1.1, v1.2, v1.3, v1.4 | frozen, tagged `prereg-v1` ... `prereg-v1.4` | `preregistration/` |
| DecodeME data, atlas, panel (19 GWAS + PGC3), FinnGen (4) | downloaded | see `data/README.md` |
| rsID map, specificity matrix, split reference | built | `scripts/build_*.py`, `split_reference.py` |
| Heritability (own LDSC), simulation tests | done; 9 tests pass | `scripts/ldsc.py`, `tests/` |
| MAGMA gene-level + models A/B for `decodeme_gwas_1` | **done, unopened** | `AAN_data/derived/magma/decodeme_gwas_1` |
| MAGMA for Alzheimer's (control) | done | `.../magma/alzheimer` |
| MAGMA for height, schizophrenia (controls) | **not done** (height reader bug, fixed in next commit; schizophrenia never started) | rerun, section 4 |
| MAGMA for `decodeme_gwas_2`, female, male, infection subsets, panel traits | not done | rerun |
| Simplified S-LDSC (own) | done for Alzheimer's, schizophrenia, RA, IBD, ME/CFS gwas_1 and gwas_2 | `results/sldsc_*.tsv` (ME/CFS ones unopened) |
| Discovery inputs (groups, HPA cell types, Reactome) | built | `AAN_data/derived/gene_covar_*.txt`, `reactome.sets` |
| Background statistics and comparison with earlier research | done | `docs/08`, `figures/fig1_background.png` |
| Report draft (methods and limits only) | draft | `paper/DRAFT_report.md` |

## 4. What to do next, in order
1. Set up the machine (`docs/09_new_machine_setup.md`) or copy `AAN_data` from the first machine.
2. Resume the heavy runs (they skip anything already finished):
   `cd scripts` then `python run_overnight.py` (ME/CFS gwas_1 is done, so it moves to the controls, the other ME/CFS GWAS,
   LDSC for all, then the panel). Run `python run_extension.py` for the FinnGen neglected conditions after it finishes.
3. When height, schizophrenia and Alzheimer's all have `cluster_B.gsa.out.txt`: `python check_controls.py`.
   Needs `G3 overall: PASS`. If FAIL, stop and debug.
4. Cell-type results: `python calibrate.py --target decodeme_gwas_1`, then `python gene_property.py --name decodeme_gwas_1 --model C`
   (after `make_conditioning_covar.py`, which needs depression, BMI and insomnia MAGMA runs).
5. Confirmation: `python confirm_robustness.py --name decodeme_gwas_1 --step locus`, then `--step permute`; S-LDSC outputs
   already exist (`results/sldsc_decodeme_gwas_1.tsv`); then `python make_evidence_ledger.py`.
6. Discovery layers: `python run_layers.py`, then `python analyze_layers.py`.
7. Translational: `python rank_targets.py --trait decodeme_gwas_1`; benchmark with `--trait migraine --disease migraine`,
   `rheumatoid_arthritis`, `ibd`; then `python smr_lite.py --candidates results/candidate_genes_decodeme_gwas_1.tsv`;
   structures with `python structure3d.py GENE1 GENE2 ...`.
8. Figures: `python make_figures.py`, 3D brain in `figures/fig7_brain3d.png` and `.html`.
9. Fill `paper/DRAFT_report.md` Results from the actual numbers, then write the abstract (under 300 words) last, in your words.
10. Application mechanics: send the three signature requests (parent or guardian, teacher, mentor) NOW if not sent; confirm the
    "no human subjects, no animal subjects" answers with your teacher; upload abstract, report, bibliography as PDF by 19 Oct.

## 5. Where data and outputs live
| Thing | Path on the first machine | Notes |
|---|---|---|
| Repo | `D:\AAN` | this repository |
| Data and big outputs | `D:\AAN_data` | about 60 GB; never committed |
| Python environment | `D:\AAN_data\venv` | recreate with `requirements.txt` |
| Reference genotypes split by chromosome | `C:\AAN_ref\by_chr` | on the SSD; made by `split_reference.py` |
| MAGMA binary | `D:\AAN_data\tools\magma\magma.exe` | v1.10 Windows build |
| All paths are hard-coded as `D:/AAN_data`, `C:/AAN_ref`, `D:/AAN` | | on another machine search and replace them in `scripts/` |

## 6. Known problems and fixes already made
- DecodeME files are raw REGENIE output; only variants in `gwas_qced.var.gz` are used (logged in `DEVIATIONS.md`).
- The standard `ldsc` package does not build on Windows (needs pysam): own implementation in `scripts/ldsc.py`.
- MAGMA is slow with the full reference; the per-chromosome reference made it about 3 times faster.
- Windows MAGMA writes `*.gsa.out.txt` and `annot.genes.annot.txt`; the scripts look for both.
- Background jobs die if the computer sleeps or the session ends. Keep the machine awake during long runs.
- The orchestrator crashed once when one trait failed; `run_panel.py` now carries on after a failed trait.
- The height GWAS file has a different column layout; `gwas_readers.read_harmonised` now handles both layouts.

## 7. Honest scores (see `docs/04_honest_grades.md`)
Advanced 5 (design), foundation 5, impact 4, novelty 3 to 4 (neuron localisation already reported by DecodeME and informally),
uniqueness 4 to 5 (provisional, re-search again on 17 Oct). Do not claim higher in the application.

## 8. Prior-art re-search still to do
Run on 17 Oct: PubMed, bioRxiv, medRxiv for "ME/CFS" with cell type, single-nucleus, LDSC, specificity, brainstem,
hypothalamus, plus whether the DecodeME paper has been published with new analyses. Update `docs/02_prior_art_and_sources.md`.
