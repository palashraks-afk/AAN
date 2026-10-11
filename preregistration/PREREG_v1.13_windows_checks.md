# Pre-registration v1.13: ten solidity checks on the six specific clusters (written 2026-10-10, before running any of them)

Target and clusters as in DISCOVERY.md A: ME/CFS decodeme_gwas_1, six clusters Amex_153, Amex_175, DLIT_152, DLIT_150, ULIT_121, Splat_402 (the pair of amygdala
clusters is the focus). Nothing here adds a cluster or changes a declared result. Each check is run once. Rules are stated now.

W1 Matched-power panel. Recompute s with the panel limited to traits that have a power comparable to ME/CFS (at most 30 clusters at FDR < 0.05 in model A:
   alzheimer, ibd, rheumatoid_arthritis, crp and every panel trait with none), and also with only the strong traits (more than 30). Rule: a cluster is "power-robust"
   if s >= 2 under the matched-power panel.
W2 False-positive rate of the T2 rule. 200 null traits made by permuting the gene rows of the cluster covariate file (all 461 columns together) within bins of
   gene size and expression, run through MAGMA model A, scored against the real panel. Report the share of null traits with at least one cluster at FDR < 0.05 and
   s >= 2, and the share of all null clusters with nominal p < 0.001 and s >= 2. The rule is acceptable if the first share is at most 5%.
W3 Jackknife over chromosomes. From the 22 leave-one-chromosome-out runs already done, compute z and s for the six clusters and a jackknife 95% interval for s.
   A cluster is "stable" if the lower bound of the interval is above 1.5.
W4 Where the cells come from. For Amex_153 and Amex_175, report the share of their cells by dissection region (from the atlas composition table), and say whether the
   majority are from amygdala dissections. Descriptive.
W5 Subsets. For each of the six clusters, the model A p-value in the infection-onset, non-infection, female and male subsets (24 tests), Bonferroni threshold 0.05/24.
   Descriptive on overlap; no claim of a difference between subsets.
W6 Amygdala region score. The region score of the amygdala for ME/CFS (mean z of its clusters weighted by cell counts) against the same score for each panel trait.
   Rule: amygdala score above the 90th percentile of the panel, and amygdala among the top 3 of 17 regions for ME/CFS.
W7 Power projection. Scale the ME/CFS cluster z by sqrt(21,620 / 15,579) (the number of cases recruited versus used, controls scaled the same way), recompute FDR over 461
   clusters for models A, B and C where the p-values are available. Report how many of the six pass FDR < 0.05 in model B and model C under the projection. A projection, not a result.
W8 Other specificity definitions. (a) rank: the percentile rank of the ME/CFS z among the panel traits' z for the same cluster; (b) trait-weighted panel: weights proportional
   to 1 minus the profile correlation with ME/CFS; (c) top decile: whether the cluster is in ME/CFS's top 10% while not in the top 10% of at least 15 of 18 panel traits.
   Rule: a cluster is "definition-robust" if it passes (a) at the 95th percentile and (b) at s >= 2.
W9 Drivers. Drop the genes with the largest contribution (z of the gene times its specificity weight in the cluster) to the Amex_153 and Amex_175 tests: top 1, top 5, top 20,
   then rerun model A. Rule: the signal does not rest on a few genes if the cluster is still nominally significant (p < 0.05) after dropping the top 20.
W10 Analysis choices. Rerun ME/CFS gene analysis and model A with SNP-to-gene windows 0/0 kb and 50/50 kb in addition to the 35/10 kb used. Report z of the six; the panel is
   unchanged. Rule: the six keep a z within 25% of their 35/10 kb value in both runs.
