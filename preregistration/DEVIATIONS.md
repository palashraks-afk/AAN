# Deviations from PREREG_v1

## 2026-10-09 - p-value floor of 1e-300 for MAGMA input (Mac rerun, before any ME/CFS cell-type result was opened)
- The analysis was rerun on a Mac. Its MAGMA build rejects p-values below about 2.2e-308 (denormal numbers) as "not a number";
  the Windows build accepted them. Height has SNPs down to 5e-324. `run_magma.py` now clips p at 1e-300 for every trait.
- Effect: only SNPs with p < 1e-300 change, and their evidence is already far beyond any threshold. No ME/CFS cluster-level
  result had been opened. Also fixed: the `.genes.out` file suffix and the PowerShell check, both Windows-only assumptions.

Each entry: date, what changed, why, and whether any ME/CFS result had been seen.

## 2026-10-08 (01:00) - C2 permutation count and S-LDSC simplification (before any result was opened)
- C2: PREREG_v1.3 asked for 5,000 gene-level permutations. Each permutation needs a MAGMA gene-property run, so 1,000
  permutations are used (smallest possible empirical p is about 0.001, enough for the p < 0.05 rule).
- C1: the standard S-LDSC needs the 75-annotation baseline model and LD scores over all reference SNPs. The first is not
  available offline here and the second would take far longer. The own implementation uses HapMap3 SNPs for LD scores and
  two control annotations (all SNPs, expressed genes +/- 100 kb). It is for confirming and ranking clusters, not for quoting
  enrichment fold-changes, and the report must say so.
- No ME/CFS cluster-level result had been opened.

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

## 2026-10-10 - code fixes after the first full run (no analysis choice changed, no extra results opened)
- Four steps crashed and were rerun after fixing: `method_benchmark.py` and `analyze_layers.py` (z was a numpy array, needs the
  cluster index), `confirm_robustness.py` (permutation columns such as `c136_p0` were dropped by the cluster-name filter), and
  `make_evidence_ledger.py` (permutation cluster ids are integers).
- `make_conditioning_covar.py`: genes MAGMA could not test for depression, BMI or insomnia had empty fields, which MAGMA rejects. They now
  get the median z of that trait (about 0). Model C was only run after this fix.
- The permutation null uses 1,000 permutations (already recorded above). Smallest possible empirical p is 0.001.
- Results were read only after the pipeline finished, and the discovery rules in DISCOVERY.md section 2 were applied as written.
