# Pre-registration v1.9: R1 re-implemented inside MAGMA (replaces the OLS version of v1.7 R1)

Written 2026-10-10 14:50, tag `prereg-v1.9`. **Why:** gate V of PREREG v1.7 failed on controls, before any ME/CFS R1 result was run:
V1 schizophrenia Spearman 0.768 (needs 0.8), and V3 as written (no neuronal cluster at FDR < 0.05 in the height residual) failed.
V3 was stricter than the main pre-registration's own height rule (G3: the three strongest clusters are non-neuronal); the height
residual's strongest clusters are fibroblasts (z about 13), which meets the G3 wording. The OLS test also disagrees with MAGMA more than
0.8 on schizophrenia, so the two tests are not interchangeable. No ME/CFS R1 or residual-profile output had been run or opened. Only
the gene family layer F1 (unchanged) and the Windows machine's T1 first look had been seen.
The failed v1.7 gate result is kept in `results/x_R1_validation.tsv` and reported as a failed first attempt.

## R1' (replaces R1 and its gate; R2 and F1 of v1.7 unchanged)
- For a target trait, the k principal components of the other panel traits' gene-level z (Kaiser rule, cap 6, computed on genes present in
  all traits, panel traits as in v1.7 R1) are added as gene-level covariates.
- MAGMA gene-property test per cluster, one-sided positive, `condition-hide = avg_all, avg_neuron, PC1..PCk` (the same form as model C).
- Specificity score s: the 461 z-values (normal quantile of p) are robustly standardised per trait; s for ME/CFS is its value minus the
  leave-one-out panel mean, divided by the panel standard deviation, as in PREREG v1 T2. The panel null uses the same procedure for
  each panel trait with the other panel traits as the PCs.
- Call "ME/CFS-specific component enriched in cluster X": FDR < 0.05 in gwas_1, p < 0.05 in the same direction in gwas_2, s >= 2.0.
- Validation gate V' (controls only, before ME/CFS): V2' Alzheimer's keeps a Microglia cluster at FDR < 0.05 under leave-one-out
  conditioning; V3' height: the three strongest clusters under leave-one-out conditioning are non-neuronal (the G3 wording).
  If V' fails, R1' is reported as a failed method and not interpreted.

## Not allowed
Changing k, the traits, the gate or the call rule after any R1' result for ME/CFS is seen.
