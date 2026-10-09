# AAN Neuroscience Research Prize 2027 — project repository

Applicant: Palash Rakshit (high school, Redlands, CA).
Prize: AAN / Child Neurology Society **Neuroscience Research Prize**, 2027 cycle.
**Deadline: Tuesday 20 October 2026, 11:59 p.m. CT** (aim to submit Mon 19 Oct).

## Working title

**Where in the nervous system does ME/CFS genetics act — and which findings point to treatment hypotheses?**
A pre-registered, specificity-calibrated cell-type map of ME/CFS (DecodeME) with a direction-aware
drug-target dossier, benchmarked on diseases whose successful drugs are known.

## Status (2026-10-09)

Analysis is running. Read **[HANDOFF.md](HANDOFF.md)** first (state, rules, exact commands), then **[TODO_PALASH.md](TODO_PALASH.md)**
(what only the applicant can do, with dates), **[SCORES.md](SCORES.md)** (honest grades) and **[DISCOVERY.md](DISCOVERY.md)**
(what could be discovered and the pre-registered rule that decides each; results not yet opened).
No ME/CFS cell-type result has been read: the control traits must pass first. Nothing in this repo is a result yet.

## Read in this order

| File | What it holds |
|---|---|
| [PLAN.md](PLAN.md) | The full plan: question, hypotheses, design, methods, timeline, decision rules |
| [preregistration/PREREG_v1.md](preregistration/PREREG_v1.md) | Frozen thresholds and tests (git tag `prereg-v1`), written before any ME/CFS result was viewed |
| [preregistration/](preregistration/) | Also: `PREREG_v1.1_extension.md` (neglected conditions), `PREREG_v1.2_discovery_layer.md` (circuits, brain vs body, pathways, comorbidity map), `PREREG_v1.3_confirmation.md` (S-LDSC, permutation null, locus drop, independent GWAS, literature convergence), `DEVIATIONS.md` (every change, dated) |
| [paper/DRAFT_report.md](paper/DRAFT_report.md) | Report draft: methods and limitations written from the code, results left empty until they exist |
| [docs/07_run_order.md](docs/07_run_order.md) | Which script runs when, and what went wrong along the way |
| [docs/01_idea_selection_log.md](docs/01_idea_selection_log.md) | Every idea considered and why each was kept or dropped |
| [docs/02_prior_art_and_sources.md](docs/02_prior_art_and_sources.md) | What already exists, with links; the gaps; what to re-search |
| [docs/03_feasibility_check.md](docs/03_feasibility_check.md) | Data files, sizes, checksums, machine, tools, what is verified vs not |
| [docs/04_honest_grades.md](docs/04_honest_grades.md) | Self-grades on the five criteria, and why no 5 on novelty/impact |
| [docs/05_aan_requirements_and_timeline.md](docs/05_aan_requirements_and_timeline.md) | Prize rules, form fields, day-by-day schedule, checklist |
| [docs/06_guardrails_and_ethics.md](docs/06_guardrails_and_ethics.md) | Patient-safety wording, human-subjects question, data-use terms |
| [data/README.md](data/README.md) | Data manifest (names, sizes, MD5) and download instructions |

## Data

No data are stored in this repository (files are 300–450 MB each). See `data/README.md`.
Data live outside the repo at `D:\AAN_data` (the C: drive has only ~13 GB free).

## Licence / reuse

Code: MIT (add LICENSE when code exists). DecodeME summary statistics are the property of the
DecodeME team; follow their terms on the OSF page.
