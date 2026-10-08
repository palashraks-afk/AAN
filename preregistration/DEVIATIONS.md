# Deviations from PREREG_v1

Each entry: date, what changed, why, and whether any ME/CFS result had been seen.

## 2026-10-07 (evening) - definition of "known targets" for the X3 benchmark
- Open Targets changed its API after the pre-registration was written (`knownDrugs` is gone, drugs now sit
  under `drugAndClinicalCandidates`). Disease ids are the current MONDO ids.
- IBD is queried through its parent node and its two children, Crohn's disease and ulcerative colitis,
  because the parent alone listed only 4 approved targets.
- A target counts as a "known winner" only if it belongs to an approved drug that has 3 or fewer
  targets. Broad drugs (for example anticonvulsants with dozens of ion-channel targets listed for
  migraine) would make the benchmark meaningless.
- Resulting lists are in `data/benchmark_targets.tsv` (migraine 24, rheumatoid arthritis 25, IBD 18
  selective targets). No gene ranking had been computed and no ME/CFS enrichment result had been seen
  when this was decided.
