# Pre-registration v1.1: extension to neglected nervous-system conditions

Written 2026-10-08, tag `prereg-v1.1`. This only adds conditions; every rule in PREREG_v1.md applies unchanged.
State of knowledge when written: the ME/CFS cell-type (gene-property) results had NOT been opened; the control traits
had not been run; only LDSC on the unfiltered ME/CFS files (discarded) and SNP-level QC checks had been seen.

## Why
ME/CFS is the anchor. To test whether the approach says anything beyond one disease, and to speak to neglected
nervous-system conditions in general, the same pipeline is applied to conditions that (a) have a public GWAS,
(b) are neglected relative to burden, and (c) had no GWAS-based cell-type analysis in my searches of 2026-10-08.

## Conditions (fixed now)
| Condition | FinnGen R13 endpoint | Cases (after QC) | Cell-type work found in search |
|---|---|---|---|
| Trigeminal neuralgia | G6_TRINEU | 2,356 | none found |
| Cluster headache (wide definition) | G6_CLUSTHEADACHE_WIDE | 1,583 | none found |
| Meniere's disease | H8_MENIERE | 3,779 | none found |
| Dystonia | G6_DYSTON | 2,142 | not searched; treat as unknown |
Dropped before running: narcolepsy with cataplexy (G6_NARCOCATA, 263 cases, far too few).
Tinnitus, restless legs and essential tremor already have published cell-type work and are not included.

## Rules
- Same gate: LDSC h2 z >= 4 or the condition is reported as excluded for low power. With 1,500 to 3,800 cases I
  expect some to be excluded, and that is a result.
- Same models A and B, same specificity score s >= 2.0 against the same panel, same FDR rule.
- Sample size for MAGMA and LDSC is set to 480,000 for all four (FinnGen R13 has about 500,000 participants; exact
  control counts are not in the files). This changes the scale of h2 but not its z-score.
- A condition is described as "no cell-type analysis found" only with the date and the search terms in
  docs/02_prior_art_and_sources.md. No "first" claims.
- Results for these conditions are exploratory next to the ME/CFS primary analysis and are labelled that way.

## Not allowed
Adding, dropping or swapping conditions after any of their results are seen.
