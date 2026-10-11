# Pre-registration v1.12: does the specificity score work? (written 2026-10-10, before running)

Why: the ME/CFS result rests on the specificity score s. So far it was checked only on three controls (G3). A reviewer will ask whether s
recovers real specificity or just ranks noise. This turns every panel trait into the target in turn, with the other panel traits as its panel.
Nothing here touches the ME/CFS result; it is a validation of the method.

V1. For each panel trait t that passes the heritability gate, compute s for all 461 clusters using the other panel traits (leave t out) and call a
    cluster "specific" if FDR < 0.05 (model A) and s >= 2 (the model B condition is not used here; it is not available for every trait).
V2. Expected-biology check, stated now: Alzheimer's (microglia clusters), rheumatoid arthritis and IBD (immune clusters: microglia, T/B cell, and other
    immune superclusters), schizophrenia (an excitatory-neuron cluster), height (non-neuronal clusters). The score "recovers" a trait if at least one specific cluster
    belongs to the expected supercluster class. Success rule: at least 4 of these 5.
V3. Base rate: for each trait, the share of FDR-significant clusters that are specific, and the median over traits. ME/CFS is then placed on that distribution
    (6 of 11 = 55%). No claim of "unusual" unless ME/CFS is above the 90th percentile of the traits.
V4. Report which cell-type classes are specific in the most traits (to show what "specific" looks like when the method works).
Whatever V1 to V4 show is reported. They cannot add or remove ME/CFS clusters.
