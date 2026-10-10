# Pre-registration v1.8: robustness of the six specific clusters (written 2026-10-10, before running anything below)

Status: **exploratory robustness checks on a result that already exists.** They can only weaken or support the six clusters declared
in DISCOVERY.md (Amex_153, Amex_175, DLIT_152, DLIT_150, ULIT_121, Splat_402). Nothing here creates a new discovery claim, and
nothing is rerun to improve a p-value.

R1. Permutation null at the originally pre-registered 5,000 permutations (v1.3 reduced it to 1,000; the seed is unchanged, so the
    first 1,000 are identical). Reported: empirical p for the top 20 clusters.
R2. Leave-one-trait-out panel: recompute the specificity score s with each of the panel traits left out in turn (and the traits
    that fail the heritability gate stay excluded). A cluster is "stable" if s >= 2 in at least 90% of the leave-one-out panels.
R3. Random half panels: 1,000 random subsets of half the panel traits. Report the share of subsets with s >= 2 for each cluster.
R4. Leave-correlated-traits-out: panel without schizophrenia, education years, BMI, neuroticism (the four most similar profiles,
    D4) and without the psychiatric/cognitive group (depression, anxiety, neuroticism, ADHD, autism, education years).

Stated before seeing the results: if a cluster is "stable" in R2 and in at least 80% of R3 subsets, it is called robust to the choice of
panel; otherwise it is reported as panel-dependent. No cluster is added to the six, whatever R2 to R4 show.
