# Pre-registration (DRAFT — freeze and `git tag prereg-v1` BEFORE any analysis)

Status: draft written 2026-10-07. Not yet frozen. Fill every `[SET]` before tagging. Values below are
proposals. Once tagged, thresholds do not move; deviations are logged in `preregistration/DEVIATIONS.md`.

## Question
Which human nervous-system cell types are enriched for ME/CFS genetic signal, and is the pattern
specific to ME/CFS?

## Data (frozen list)
ME/CFS: DecodeME gwas_1, gwas_2, gwas_1_female, gwas_1_infectious_onset, gwas_1_non_infectious_onset
(gwas_1_male gated by h2 rule). Atlas: Human Brain Cell Atlas v1.0 aggregated loom. Partition level:
**[SET: cluster (461) as primary; supercluster as secondary]**. Comparison-trait panel: **[SET list of
≥ 20 traits with source, N, URL before any enrichment is run]**.

## Gates (decided before looking at enrichment)
- G1 heritability: analyse a GWAS only if LDSC h2 z ≥ **[SET 4]**.
- G2 variant QC: INFO ≥ **[SET 0.8]**, MAF ≥ **[SET 0.01]**, MHC excluded.
- G3 positive controls must behave: schizophrenia enriched in cortical neurons and Alzheimer's in
  microglia at FDR < 0.05, and height NOT enriched in neural clusters. If any fails, stop and debug
  the pipeline; do not interpret ME/CFS.

## Primary tests
- T1 Enrichment per cluster: MAGMA gene-property (continuous specificity), controlling for mean
  expression; BH-FDR across clusters; significant if FDR < **[SET 0.05]**.
- T2 Specificity: a cluster is "ME/CFS-specific" iff ALL of: (a) T1 significant; (b) its ME/CFS
  enrichment statistic exceeds the **[SET 95th]** percentile of the same cluster's statistic across
  the panel traits; (c) enrichment remains significant (p < **[SET 0.05]** after FDR) conditional on the
  generic neuronal factor (mean expression across neurons).
- T3 Conditioning: repeat T1–T2 conditioning on gene-level association of depression, BMI and insomnia
  (MAGMA conditional model). Report whether each specific cluster survives.
- T4 Internal replication: T1–T2 repeated on gwas_2. Note shared cases; call it "robustness", not
  independent replication.
- T5 Subtype/sex: difference in enrichment effect size (infection vs non-infection; female vs all)
  with bootstrap 95% CI **[SET 1,000 resamples of genes/blocks]**. Descriptive.

## Translational tests
- X1 Candidate genes: top **[SET 10]** genes per specific cluster by (MAGMA Z × specificity).
- X2 Evidence tiers: A = colocalisation PP.H4 ≥ **[SET 0.8]** + MR direction consistent across
  instruments + existing drug with matching mechanism direction + LOEUF not in lowest decile;
  B = three of those four; C = candidate only.
- X3 Benchmark: on **[SET ≥ 3 positive-control diseases with known drug targets, e.g. migraine,
  rheumatoid arthritis, multiple sclerosis]**, the pipeline must place known approved-drug targets at
  a median rank percentile better than **[SET 80th]** versus size/length-matched random gene sets
  (permutation p < 0.05). If not met, all ME/CFS tiers are labelled "exploratory, uninterpreted".

## Outcomes that will be reported whatever happens
- Null: no ME/CFS-specific cluster → "ME/CFS signal is indistinguishable from the generic neuronal
  signal of polygenic brain traits at this sample size."
- Excluded analyses with the gate that excluded them.
- Sensitivity: top-5% / 10% / 20% specificity gene sets; panel leave-one-out; atlas partition level.

## What would be a deviation
Changing a threshold, adding/removing a panel trait after seeing ME/CFS results, switching the
partition level, or re-defining "specific". Each logged with date and reason.

## Stopping rule
If the go/no-go gate on 13 Oct is failed, report T1–T2 + limits only and drop X-tests and stretch.
