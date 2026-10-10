# Pre-registration v1.10: independent expression source (GTEx bulk tissues) for the amygdala finding, plus peripheral nerve and pituitary

Written 2026-10-10 (Mac extension pipeline), tag `prereg-v1.10`. Adds analyses only.
State of knowledge when written: the main pipeline declared amygdala excitatory neurons (Amex_153, Amex_175) as the only panel-robust specific
clusters; the whole-body HPA layer showed somatotrophs and lactotrophs at FDR < 0.05 (pituitary cell types), and the Siletti atlas has no
peripheral nerve. No GTEx-based ME/CFS result has been computed. The three hypotheses below come from those existing results and the atlas
gap, and are fixed now.

## Why
The amygdala result rests on one atlas (three donors, single-nucleus). GTEx is a different dataset (bulk RNA-seq, other donors, other
technology) that includes a Brain-Amygdala sample. Agreement would be an independent check; disagreement would weaken the claim. Nerve - Tibial
is a peripheral nerve, so it tests part of the gap the brain atlas cannot cover. Pituitary follows from the HPA somatotroph/lactotroph signal.

## Data
GTEx v8 gene median TPM by tissue (54 tissues, 13 of them brain; `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`).
Genes mapped to Entrez ids by symbol using NCBI38.gene.loc; genes with mean TPM >= 1 across tissues kept.

## Specificity scores
- Across-tissue (S54): per gene, TPM in the tissue divided by the sum over the 54 tissues, rank-normalised within the tissue column.
- Within-brain (SB): per gene, TPM in a brain tissue divided by the sum over the 13 brain tissues, rank-normalised.
Covariates: avg54 (log10 mean TPM over all 54) and avgbrain (log10 mean over the 13 brain tissues).

## Tests (MAGMA gene-property, one-sided positive, on genes.genes.raw)
- H1 Amygdala: SB for Brain - Amygdala, conditioning on avg54 and avgbrain (asks whether amygdala-specific genes carry more ME/CFS signal
  than the genes specific to other brain regions).
- H2 Peripheral nerve: S54 for Nerve - Tibial, conditioning on avg54.
- H3 Pituitary: S54 for Pituitary, conditioning on avg54.
- Exploratory: all 54 tissues with S54 (BH-FDR across 54), and all 13 brain tissues with SB (BH across 13).
- Targets: decodeme_gwas_1 (primary), decodeme_gwas_2 (same cases, different controls), and the 19 panel traits.

## Calls
- H1/H2/H3 "supported": p < 0.05/3 in gwas_1, p < 0.05 in the same direction in gwas_2, and a specificity score s >= 2 against the panel
  (each trait's z across the tested tissues robustly standardised; s = ME/CFS minus the panel mean for that tissue, over the panel SD).
- Anything from the exploratory sets is "observed, not supported" unless it passes FDR < 0.05 in gwas_1 and p < 0.05 in gwas_2.
- A null on H1 is reported as: the amygdala result does not replicate in an independent expression source.

## Not allowed
Adding tissues or hypotheses, changing covariates or thresholds after any ME/CFS GTEx result is seen.
