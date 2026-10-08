# Pre-registration v1 (frozen 2026-10-07, before any ME/CFS enrichment result was viewed)

Tag: `prereg-v1`. Anything that changes after this tag goes in `DEVIATIONS.md` with a date and a reason.
State of knowledge at freezing: the DecodeME summary statistics are downloaded, the rsID map and the
specificity matrix exist, MAGMA for GWAS-1 was running and its output had not been opened. LDSC had
only been run on schizophrenia (PGC3) as a code sanity check.

## Question
Which human nervous-system cell types carry ME/CFS genetic signal, and is any of that signal specific
to ME/CFS rather than shared by polygenic brain traits in general?

## Data (fixed)
- ME/CFS: DecodeME `gwas_1` (primary), `gwas_2` (robustness; shares all cases with gwas_1),
  `gwas_1_female`, `gwas_1_infectious_onset`, `gwas_1_non_infectious_onset`, `gwas_1_male`.
- Atlas: Human Brain Cell Atlas v1.0, aggregated loom, **461 clusters** (primary level), 31
  superclusters (secondary level).
- Gene sets: 14,049 genes expressed (mean CPM >= 1 across clusters) and present in the NCBI38 gene file.
- Comparison panel (the null for specificity): the 23 traits in `data/panel.tsv` that are downloaded and
  pass the heritability gate, plus PGC3 schizophrenia. Chronic pain, fibromyalgia, type 1 diabetes
  are excluded because no suitable full-genome file was found; this is decided before looking at results.

## Gates
- **G1 heritability.** A GWAS is analysed only if LDSC h2 z-score >= 4 (own implementation, tested by
  simulation). Males will probably fail and are then reported as excluded.
- **G2 variants.** MAF >= 0.01; autosomes only; a SNP must match a 1000G EUR rsID by GRCh37 position and
  allele pair. No separate INFO filter (the DecodeME release is already QC'd; to be confirmed against
  `gwas_qced.var.gz`).
- **G3 positive/negative controls (operational definitions).** Pipeline is trusted only if all hold:
  (i) for PGC3 schizophrenia at least one cluster in an excitatory-neuron supercluster has FDR < 0.05
  in model A; (ii) for Alzheimer's (Bellenguez) at least one Microglia cluster has FDR < 0.05;
  (iii) for height (Yengo) the three highest-z clusters are all non-neuronal superclusters.
  If G3 fails, stop and debug; ME/CFS results are not interpreted until it passes.

## Models (MAGMA 1.10, gene-property, one-sided positive)
- Model A: specificity (rank-normalised) of each cluster, conditioning on log10 mean expression across
  all clusters (`avg_all`).
- Model B: model A plus conditioning on log10 mean expression across neuronal clusters (`avg_neuron`).
- Model C (T3): model B plus the gene-level z-scores of depression (broad), BMI and insomnia GWAS.

## Primary tests
- **T1 localisation.** ME/CFS gwas_1, model A, Benjamini-Hochberg across the 461 clusters, FDR < 0.05.
- **T2 ME/CFS-specific cluster.** A cluster is called ME/CFS-specific only if ALL hold:
  (a) T1 significant; (b) model B p < 0.05 for that cluster; (c) specificity score s >= 2.0, where for each
  trait the 461 cluster z-values (z = normal quantile of the one-sided p) are robustly standardised
  (subtract the median, divide by 1.4826 x MAD across clusters), and s is ME/CFS's standardised value
  minus the panel mean for that cluster, divided by the panel standard deviation.
  Panel percentile is reported but is not part of the rule.
- **T3 conditioning.** For each T2 cluster: does model C p stay < 0.05?
- **T4 robustness.** The same T1/T2 on gwas_2.
- **T5 subtype and sex (descriptive).** Standardised-profile differences between infectious-onset and
  non-infectious-onset GWAS, and between female-only and gwas_1, reported with the nominal MAGMA
  beta/SE difference z. Samples overlap (shared controls and cases) so these carry no inferential claim.

## Translational layer
- **X1 candidates.** For each T2 cluster, the 10 genes with the highest value of
  (MAGMA gene z) x (cluster specificity).
- **X2 evidence tiers.** A: colocalisation PP.H4 >= 0.8 with a brain eQTL/pQTL, MR direction consistent
  across instruments, an existing drug whose mechanism direction matches, and LOEUF not in the lowest
  decile. B: three of the four. C: candidate only.
- **X3 benchmark (must pass or the ME/CFS tiers are labelled exploratory).** Run the same ranking on
  migraine, rheumatoid arthritis and inflammatory bowel disease using the panel GWAS. Known targets of
  approved drugs (Open Targets, approved-drug status read before ranking) must have a median rank
  percentile >= 80 and permutation p < 0.05 against matched random gene sets (same size, matched on
  gene length and expression). Fewer than 3 diseases passing = exploratory.

## Reported whatever happens
A null T2 (no specific cluster) is a result. Excluded GWAS and the gate that excluded them. Sensitivity:
top-5% and top-20% specificity sets; leave-one-out over panel traits; panel without the related-condition
traits; supercluster level.

## Not allowed after freezing
Changing a threshold, panel membership, cluster level, or the definition of "specific" after seeing ME/CFS
results.

## Stopping rule
If the core T1/T2 pipeline is not producing results by the end of 13 Oct, drop the translational layer and
the stretch panel.
