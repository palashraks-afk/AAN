# Pre-registration v1.6: related conditions outside the null panel

Written 2026-10-08, tag `prereg-v1.6`. No ME/CFS cluster-level result opened.

## Added
- Long COVID, strict case definition (Long COVID after test-verified infection) against population controls, global
  multi-ancestry meta-analysis of 12 GWAS (Lammi et al., Nature Genetics 2024; GWAS Catalog GCST90454540, harmonised, GRCh38).
- UK Biobank-coded chronic fatigue syndrome (GCST90038694), already in the plan as a low-power independent check (v1.3 C4).

## How they are used
- Same heritability gate (LDSC z >= 4). If a file fails the gate it is reported as excluded for low power.
- If it passes: MAGMA models A and B and S-LDSC like every other trait, but it is NEVER part of the null comparison panel
  (`calibrate.NOT_IN_PANEL`).
- Used in the cellular comorbidity map (D4, v1.2) and the profile-similarity figure, and in concordance with ME/CFS (C4).
- Caveat to state: the long COVID file is a multi-ancestry meta-analysis analysed against a European LD reference, so its
  heritability and cell-type results are approximate.

## Not available, so not used
The 2026 fibromyalgia GWAS (Nature Communications 2026; 23,209 European cases) lists no full summary statistics in the GWAS
Catalog yet. The larger Nature Medicine fibromyalgia paper promises release on publication. Re-check on 17 Oct; if they
appear they will be added under this same rule (related condition, outside the null).

## Not allowed
Putting any of these into the null panel, or adding conditions after seeing ME/CFS results.
