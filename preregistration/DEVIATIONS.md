# Deviations from PREREG_v1

Each entry: date, what changed, why, and whether any ME/CFS result had been seen.

## 2026-10-08 (00:00) - G2: the DecodeME files are NOT already QC'd
- PREREG_v1 G2 assumed the shared GWAS files were already QC'd. They are not: the README
  (`README_shared_sum_stats.txt`) says they are raw REGENIE output filtered on MAF only, and tells users to
  keep only the variants in `gwas_qced.var.gz` (MAF >= 0.01, INFO >= 0.9 and a DENTIST-style test).
- Found because the first MAGMA gene table had implausible top genes (p of 1e-47) driven by hundreds of
  SNPs in the pericentromeric regions of chr15 (22.0-22.9 Mb) and chr21 (9.04 Mb), far stronger than
  DecodeME's own published loci.
- Change: every DecodeME analysis now keeps only variants on the QC list. The first LDSC and MAGMA runs
  (without the filter) are discarded and redone.
- What had been seen: the 15 most significant SNPs and the top genes of the unfiltered run, as a data
  check. No cluster-level (cell-type) result had been opened.

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
