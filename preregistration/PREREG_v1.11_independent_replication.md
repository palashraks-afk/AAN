# Pre-registration v1.11: independent single-nucleus atlas and an independent fatigue cohort

Written 2026-10-10 (Mac extension pipeline), tag `prereg-v1.11`. Adds analyses only.
State of knowledge: the main pipeline's six specific clusters (Amex_153, Amex_175, DLIT_152, DLIT_150, ULIT_121, Splat_402) and the GTEx non-replication
(v1.10) are known. No ME/CFS result has been computed on the data below. I looked at the cell-type labels and counts of the atlas before writing this
(labels only), and at nothing from ME/CFS.

## Why
The amygdala result rests on one atlas and one GWAS family (DecodeME). Two independent checks: (1) a different single-nucleus dataset with amygdala and
striatal cell types (Tran et al. 2021, LIBD, CELLxGENE collection 22ccff12-9a18-4b01-ba10-302df57a01f8; 5 donors for amygdala, 8 for nucleus accumbens);
(2) an independent cohort with no shared participants, the FinnGen R13 endpoint R18_MALAI_FATIG (malaise and fatigue; 31,709 cases).

## P1 Independent atlas (Tran 2021)
- Data: the AmyG and NAc datasets only (the other three regions were not downloaded; fixed). Raw counts aggregated per `author_cell_type`; groups with fewer
  than 50 nuclei and groups named `drop.*` (doublets, low quality) excluded. Result: 35 groups (15 AmyG, 20 NAc). Specificity as in v1.10 (fraction of
  expression across groups, rank-normalised), conditioning on mean expression, MAGMA gene-property, one-sided positive.
- Primary hypotheses (a priori, from the main pipeline): H-A1: the amygdala excitatory groups AmyG Excit_A (344 nuclei) and AmyG Excit_C (55 nuclei).
  Supported if, for either, p < 0.025 (Bonferroni over 2) in gwas_1, p < 0.05 in gwas_2 in the same direction, and specificity score s >= 2 against the 19-trait panel.
- Secondary (exploratory, no confirmation language): H-A2 the medium spiny neuron groups of NAc (D1 and D2) as the eccentric-MSN story of earlier reports;
  all 35 groups with BH-FDR.
- Stated limitations: very few amygdala excitatory nuclei (443 glutamatergic in the dataset; Excit_B has 44 and is excluded), coarse labels, 5 donors. A null here is
  weak evidence against the amygdala result and is reported as uninformative about it, not as a refutation.

## R1 Independent cohort (FinnGen R18_MALAI_FATIG)
- Gate: LDSC heritability z >= 4; otherwise reported as excluded for low power.
- MAGMA and the Siletti 461-cluster model A as for every other trait, with n = 480,000 as in v1.1.
- H-R: the mean z of the six specific clusters in the FinnGen fatigue profile exceeds the 95th percentile of the mean z of 10,000 random sets of six clusters
  drawn from the same profile. Also reported: how many of the six have z > 0, and the 11 main clusters' z.
- Caveat stated now: "malaise and fatigue" is a symptom code, broader than ME/CFS, so a null does not refute the ME/CFS result and a positive is support, not replication.

## Not allowed
Adding groups, regions or cohorts, or changing thresholds, after any ME/CFS result on these data is seen. A fetal/peripheral atlas test (Cao et al. 2020) will be
pre-registered separately (v1.12) after its labels are inspected.
