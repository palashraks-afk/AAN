# Pre-registration v1.2: discovery layer on top of the known neuronal signal

Written 2026-10-08, tag `prereg-v1.2`. Adds exploratory analyses only; PREREG_v1 and v1.1 stay as they are.
State of knowledge when written: ME/CFS gene-level results seen as a QC check (top genes match the published loci);
no cluster-level ME/CFS result opened; control traits still running.

## Why
DecodeME's team and others already say the ME/CFS signal sits in neurons. That is the starting point. These layers ask
what the known result does not: which kind of neuron and transmitter system, whether the signal is brain-specific or
body-wide, which biological programmes sit behind it, and which other conditions share the same cellular profile.
All four reuse the finished SNP-to-gene step (MAGMA `genes.genes.raw`) so they are cheap, and all are run on the
same comparison panel so each result can be calibrated the same way.

## D1 circuit resolution inside neurons
- Groups from the atlas annotation: transmitter class (GABAergic, glutamatergic, GABA-glutamate mixed, glycinergic,
  monoaminergic or cholinergic, unannotated) and the first token of the subtype annotation (for example MSN-D1,
  MSN-D2, INT-SST, INT-PVALB, INT-VIP, DG-GRAN). A group needs at least 3 clusters.
- Group specificity: cell-weighted mean CPM over the group's clusters, divided by the sum over groups, rank-normalised.
- Test: MAGMA gene-property, conditioning on mean expression (model A form) and on mean neuronal expression (model B form).

## D2 brain versus body
- Human Protein Atlas single-cell-type consensus expression (all organs). Same specificity construction, same tests.
- Reported: rank of every cell type, and whether brain neurons beat immune, muscle and other cell types.
- A claim of "brain-specific" requires the best neuronal cell type to beat the best non-neural cell type at FDR < 0.05 in
  both gwas_1 and gwas_2.

## D3 biological programmes
- Reactome pathways (10 to 500 genes). MAGMA competitive gene-set test, conditioning on mean expression, on mean
  neuronal expression, and in a second model on the specificity of the single ME/CFS cluster with the smallest model-A p.
- FDR across pathways. A pathway "adds beyond the cell type" only if it stays FDR < 0.05 in the second model.

## D4 cellular comorbidity map
- Spearman correlation of the 461-cluster z-profile between ME/CFS and every panel trait, and between every pair of
  traits, then hierarchical clustering. Reported as a map of which conditions share ME/CFS's cellular profile.
- Null: profiles of ME/CFS built from permuted gene-level z within matched gene-size bins (1,000 permutations) to say
  which correlations are larger than chance.

## Rules for calling anything "supported"
- BH-FDR within each layer, never across layers.
- A finding is "supported" only if it also reaches p < 0.05 in the same direction in `decodeme_gwas_2` (robustness; shares
  cases) and, for D1 and D2, ME/CFS exceeds the panel (specificity score s >= 2.0 as in PREREG_v1).
- Everything else is reported as "observed, not supported".
- The control gate G3 must pass before any of this is read.

## Not allowed
New groups, new pathway collections, new thresholds, or dropping a layer after results are seen.
