# Self-review log

Problems found by checking my own work, what fixed them, and what could still be wrong. Updated 2026-10-08.

## Found and fixed
| # | Problem | How it was found | Fix |
|---|---|---|---|
| 1 | Pre-registration assumed the DecodeME files were already quality-controlled | implausible top genes (p of 1e-47 in pericentromeric chr15 and chr21) | filter to `gwas_qced.var.gz`; logged in DEVIATIONS.md; first runs discarded |
| 2 | `--merge` of MAGMA batches failed with absolute paths | non-zero exit status | run it from inside the folder with relative names; stitch the per-chromosome gene tables |
| 3 | MAGMA took hours | timing test on one chromosome | per-chromosome reference files on the SSD, `--batch N chr` |
| 4 | Panel readers crashed on 7 of 20 GWAS files (a second column layout) | height run failed | `read_harmonised` handles both layouts, refuses files with too few usable SNPs |
| 5 | One failed trait stopped the whole panel run | orchestrator log | `run_panel.py` logs the failure and continues |
| 6 | The null comparison panel would have included the neglected conditions and the UK Biobank fatigue GWAS | code review before any ME/CFS result was read | `calibrate.NOT_IN_PANEL`, used in calibration, layers and figures |
| 7 | Output file names from the Windows build of MAGMA (`.gsa.out.txt`) did not match what the loaders expected | first end-to-end test | loaders try both names |
| 8 | OSF and the GTEx, FinnGen, Open Targets, GWAS Catalog interfaces had changed or rate-limited | download failures | retry with backoff; new API queries; documented in docs/07 |
| 9 | Background jobs died when the session ended or the computer slept | process list showed nothing running after 14 hours | detached launch; scripts skip finished work; keep the machine awake |

## Checks that passed
- Own LDSC: recovers simulated heritability and intercept (4 tests); gives an ME/CFS liability-scale value close to DecodeME's published 0.095 at a prevalence of 0.5%.
- Own S-LDSC: unit tests pass (5); recovers expected cell classes for schizophrenia (excitatory neurons), Alzheimer's (microglia),
  rheumatoid arthritis and IBD (immune cells). See `results/controls_sldsc.md`.
- ME/CFS gene-level results after the QC filter match DecodeME's published loci (ARFGEF2, CSE1L, STAU1, TAOK3, histone cluster).

## Still to check or could be wrong
- VERIFIED 2026-10-08 17:45: the per-chromosome reference gives identical gene-level results to the full reference on
  chromosome 22 of decodeme_gwas_1 (444 genes in both, maximum |z difference| 0, correlation 1.000000, same SNP counts;
  `scripts/check_split_reference.py`).
- MAGMA-based controls for height and schizophrenia have not finished; nothing about ME/CFS is read until they pass.
- The 1000G reference lacks the synonym file in per-chromosome form, so panel GWAS using merged rsIDs lose a few SNPs.
- S-LDSC uses HapMap3 SNPs only and two control annotations, not the 75-annotation baseline: good for confirming and ranking,
  not for quoting enrichment sizes.
- Direction of effect from SMR without HEIDI cannot separate causality from linkage.
- Panel GWAS differ in size and design; many use UK Biobank, so controls overlap with DecodeME's controls.
- Case-control sample size for FinnGen conditions was set to a round 480,000 (affects the scale of h2, not its z-score).
- The literature convergence region list was fixed from searches; other reviewers might choose a slightly different list.
- Prior-art search is time-limited. Re-run it on 17 Oct.
