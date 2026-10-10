# Pre-registration v1.7: the ME/CFS-specific genetic component, and gene families

Written 2026-10-10 on the Mac extension pipeline, tag `prereg-v1.7`. Adds analyses only; no earlier rule changes.
State of knowledge when written: control gate G3 passed on the Windows machine and a first look at ME/CFS localisation (T1, model A,
11 clusters, none in model B) is recorded in DISCOVERY.md. **On the Mac, no ME/CFS cluster-level table was opened.** Only
gene-level files exist for ME/CFS here. Nothing below was designed from the shape of any ME/CFS cluster result.

## Why
PREREG v1 asks whether a cluster's ME/CFS signal is bigger than other traits show for that cluster (the s score). That is a
test at the cluster level. It cannot say how much of ME/CFS's gene-level signal is left once the genetics that all the
other traits share is removed. R1 asks that directly at the gene level and then asks where the leftover signal sits. F1
asks which classes of drug-relevant genes carry signal.

## R1 The ME/CFS-specific component (gene-level residualisation)
Inputs: MAGMA gene-level z (ZSTAT) for ME/CFS `decodeme_gwas_1` and for every panel trait that passes G1 and was run on the same
gene set (not ME/CFS, not the related-condition or neglected-condition sets of v1.1 and v1.6).
1. Matrix Z of genes x panel traits (genes present in all). Standardise each column. PCA on the correlation matrix. Keep
   k = number of PCs with eigenvalue >= 1, capped at 6 (fixed rule; sensitivity at k = 3 and k = 6 only).
2. Regress ME/CFS gene z on the k PCs, log(NSNPS) and log(NPARAM). The residual is the ME/CFS-specific gene signal r.
3. Cell-type test: for each of the 461 clusters regress r on the cluster's rank-normalised specificity, with `avg_all` and
   `avg_neuron` as covariates (OLS, one-sided positive). Because OLS ignores gene-gene LD, the p-value is taken from a
   permutation of r within 10 x 10 bins of gene size and mean expression (2,000 permutations; empirical p, add-one).
   BH-FDR across the 461 clusters.
4. Specificity score: for each trait the 461 residual z-values are robustly standardised (median, 1.4826 x MAD); s for ME/CFS is
   its value minus the leave-one-out panel mean for that cluster, divided by the panel standard deviation, as in PREREG v1 T2.
5. Panel null: the same procedure with each panel trait taken in turn as the target and the other panel traits as the PCs.

### Validation gate V (decided before running; if it fails, R1 is reported as a failed method and not interpreted)
- V1: on schizophrenia, Alzheimer's and height, the OLS-plus-permutation cluster z-profile (raw gene z, no residualisation)
  correlates with the MAGMA gene-property z-profile (model A) with Spearman >= 0.8.
- V2: the residual of Alzheimer's (leave-one-out) still has a Microglia cluster at FDR < 0.05, because microglial biology is
  specific to it in the panel. If residualisation erases it, the method removes real specific signal.
- V3: height residual has no neuronal cluster at FDR < 0.05 (negative control).

### Calls
- "ME/CFS-specific component enriched in cluster X": FDR < 0.05 in R1 for gwas_1, p < 0.05 in the same direction in gwas_2
  (the same residualisation applied to gwas_2), and s >= 2.0.
- "Nothing specific remains": zero such clusters. This is a reportable result: it says the cellular signal is shared with the panel.
- Also reported (descriptive): how much of ME/CFS gene z variance the k PCs explain, and the number of ME/CFS-significant
  clusters in the raw profile that lose significance after residualisation.

## R2 Genes carrying the specific component
- Genes with residual z >= 3.5 (about p < 2e-4), their overlap with the DecodeME loci, and their Open Targets approved-drug
  status. Hypotheses for researchers only; the benchmark rule of PREREG v1 X3 applies before any target is interpreted.

## F1 Gene families
- Families fixed now, by HGNC gene-group name matching (case-insensitive): "channel", "G protein-coupled receptor",
  "neuropeptide|peptide hormone", "solute carrier", "kinase", "nuclear receptor", "cadherin", "immunoglobulin". A family needs
  at least 30 genes.
- MAGMA competitive set test on `genes.genes.raw` for gwas_1 and gwas_2, conditioning on `avg_all` and `avg_neuron`.
  BH-FDR across the 8 families. "Supported": FDR < 0.05 in gwas_1 and p < 0.05 in gwas_2, same direction.

## Reporting
Evidence ledger entries for R1, R2 and F1 with the numbers, the validation gate result and the panel-null values. A failed gate
is reported as such. No "first" or "never" wording; the 17 Oct prior-art search decides what can be said about novelty.

## Not allowed
Changing k, the family list, the thresholds, the traits in the panel or the validation gate after any R1 or F1 result is seen.
