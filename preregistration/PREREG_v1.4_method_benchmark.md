# Pre-registration v1.4: benchmark of the cell-type methods on control traits

Written 2026-10-08, tag `prereg-v1.4`. Descriptive method comparison; it does not change any earlier rule.
State of knowledge: no ME/CFS cluster-level result opened. Controls were still running.

## Methods compared
1. MAGMA gene-property with continuous specificity (model A of PREREG_v1).
2. MAGMA competitive set test on the top 10% most specific genes of each cluster (`top10_sets.txt`), conditioning on
   gene size and the gene-level covariates MAGMA uses by default.
3. Simplified S-LDSC (`sldsc_lite.py`, PREREG_v1.3 C1).

## Traits
Schizophrenia (expect excitatory cortical neurons), Alzheimer's disease (microglia), rheumatoid arthritis and type 1
diabetes if available (immune clusters), height (non-neural, expect no neuronal enrichment), plus the rest of the panel
for agreement statistics. ME/CFS is NOT part of this benchmark.

## Metrics, fixed now
- Recall: for each positive control, is at least one expected cluster FDR < 0.05 under each method (yes/no).
- False positives: number of neuronal clusters with FDR < 0.05 for height under each method.
- Agreement: Spearman correlation of the 461-cluster z-profile between each pair of methods, per trait; median over the panel.
- Cost: wall-clock seconds per trait.
- Expected clusters are defined from the atlas supercluster labels: excitatory cortical neurons (upper-layer and deep-layer
  intratelencephalic, near-projecting, corticothalamic and 6b), Microglia, and for immune traits the Miscellaneous
  supercluster (B, T, NK, monocyte clusters).

## Reporting
A table of the metrics, no ranking of methods beyond the table, and a plain statement of any method that fails a control.
If a method fails the schizophrenia or Alzheimer's control it is not used for confirmation of ME/CFS.
