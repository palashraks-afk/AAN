# Pre-registration v1.14: meta-analysis cluster profile, replicated-loci genes, and related conditions (tasks M1, M3, M4 of NEXT_20)

Written 2026-10-10 (Mac session), tag `prereg-v1.14`. Adds analyses only; no earlier rule changes. Nothing below has been computed on ME/CFS or on the related conditions.

## M1 DecodeME + MVP meta-analysis (descriptive, NOT independent)
Source: Maccallini et al. 2026, Zenodo record 20204356 (CC-BY 4.0), `GWAS_METAL_DME_1_MVP_GRCh38.tsv.gz`: a meta-analysis that contains the DecodeME GWAS-1 cases, so it cannot test replication.
Run MAGMA and the Siletti 461-cluster model A exactly as for gwas_1. Report the z and p of the six specific clusters and of the 11 main clusters next to their DecodeME-only values.
Stated reading rule: if the six clusters' mean z rises, adding 3,891 MVP cases does not dilute the signal; if it falls, the extra cases add noise or a different phenotype. Their paper's own Table 5 shows
the MVP cohort alone has near-zero tissue effects, so a fall is expected. No claim of replication is allowed from this layer.

## M3 Replicated-loci genes
Gene set fixed now: CLYBL, BICD1, GRIN2A, CSMD1, RORA (named for the biobank study of Slaughter et al. 2026; the other replicated genes are not available to me and are not added).
Statistic: the mean rank-normalised specificity of these genes in each of Amex_153 and Amex_175, against 10,000 random sets of five genes matched on the decile of mean expression.
"Supported" if the empirical one-sided p < 0.025 (Bonferroni over the two clusters) for either. Otherwise "not supported". Five genes is a very small set; a null is uninformative.

## M4 Related conditions: is the amygdala specificity also present?
Three conditions, fixed now: irritable bowel syndrome (EBI harmonised file already in the panel, run as a target), FinnGen R13 M13_FIBROMYALGIA (4,005 cases), FinnGen R13 PAIN
(246,393 cases, "Pain (limb, back, neck, head abdominally)"). (Long COVID has no open full-genome file that I could obtain in time and is not tested; malaise and fatigue was reported under v1.11.)
Rule as in v1.11 H-R: the mean z of the six specific clusters (Amex_153, Amex_175, DLIT_152, DLIT_150, ULIT_121, Splat_402) in the condition's 461-cluster profile exceeds the
95th percentile of 10,000 random sets of six clusters from the same profile; Bonferroni over the three conditions (p < 0.05/3 = 0.0167). Heritability gate z >= 4 first; a failing condition is reported as excluded.
These conditions never enter the null panel. Reading: a positive means the six clusters are shared with related conditions (support for a general fatigue/pain/neuronal biology), not specific to ME/CFS.

## Not allowed
Adding genes, conditions or clusters, or changing thresholds, after any result of these layers is seen.
