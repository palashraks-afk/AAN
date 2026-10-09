# Pre-registration v1.5: independent signals, region-level tests, 3D landscape, per-locus contributions

Written 2026-10-08 (evening), tag `prereg-v1.5`. Adds exploratory layers. Nothing earlier changes.
State of knowledge: no ME/CFS cluster-level result opened. Alzheimer's control passes rule (ii) of PREREG_v1 G3 in MAGMA
(microglia, FDR 3e-7); schizophrenia and height MAGMA runs were still going.

## L1 How many independent cell-type signals? (forward selection)
- Start from the cluster with the smallest model-A p-value. Add it as a conditioning variable (MAGMA `condition`), retest all
  other clusters, take the best remaining one if FDR < 0.05, and repeat, up to 5 independent clusters.
- Reported: the ordered list, each cluster's z before and after conditioning, and the superclusters involved.
- Interpretation rule: if after the first cluster nothing else survives, the signal is "one cell-type axis" (the many
  enriched clusters are redundant because they share genes); if more survive, report them as separate axes.
- The same procedure is run on the three positive controls (schizophrenia, Alzheimer's, rheumatoid arthritis) to show how many
  independent axes a well-understood trait has.

## L2 Region-level test (direct, not via cluster z-values)
- Pool cluster expression into the 17 atlas regions (cells from each dissection mapped to its region), build gene specificity
  across regions (same construction as clusters), and run MAGMA gene-property, conditioning on mean expression.
- Calibrated against the panel exactly as for clusters (specificity score s, and replication in gwas_2).
- This gives region statistics with p-values; the 3D brain figure uses these region z-values, not the cell-weighted averages.

## L3 3D cellular landscape (visual, no test)
- The 461 clusters are placed in 3D by principal components of log expression. Colour is the ME/CFS enrichment z, size is
  log cell count, hover text gives name, supercluster and z. A second version is coloured by the specificity score s.
  Purely descriptive: it shows whether enriched clusters are neighbours (one cell family) or scattered.

## L4 Per-locus contributions
- For the three clusters with the smallest model-A p, drop each genome-wide significant locus in turn (genes within 1 Mb of
  the lead variant, loci merged when closer than 1 Mb) and report the change in z. A cluster is "single-locus-driven" if
  removing one locus drops its z below the FDR threshold.

## L5 Central versus peripheral nervous system context (reporting rule)
- The whole-body layer (D2, PREREG v1.2) reports explicitly: brain neurons, glia, Schwann cells (peripheral nerve), enteric
  cells, retinal neurons and neuroendocrine cells, so the report can say whether the signal is restricted to the brain.

## Rules
- BH-FDR within each layer. "Supported" needs gwas_2 replication in the same direction at p < 0.05 (as in v1.2).
- L1 to L4 are exploratory and labelled so; none enters the T2 definition.

## Not allowed
Changing the number of axes, the conditioning order rule, the region mapping, or the loci rule after results are seen.
