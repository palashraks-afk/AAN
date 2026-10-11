# Pre-registration v1.18: is the six-cluster pattern still ME/CFS-specific once pain traits join the comparison panel, and do the two machines agree?

Written 2026-10-10 (Mac session), tag `prereg-v1.18`, before the computation. Nothing here changes the earlier headline; it quantifies a limitation (no pain trait in the comparison panel) and a reproducibility question.

## P1 Pain-augmented panel (sensitivity analysis)
Question: v1.15 M4 showed the six specific clusters (Amex_153, Amex_175, DLIT_152, DLIT_150, ULIT_121, Splat_402) carry signal in FinnGen pain. The original panel (19 traits) contains no pain trait.
Method: recompute the specificity score s for all 461 clusters exactly as in PREREG_v1 T2 (each trait's 461 z-values from MAGMA model A, robustly standardised; s = ME/CFS value minus the panel mean, divided by the
panel SD), once with the original 19-trait panel and once with the panel plus FinnGen pain and FinnGen fibromyalgia (21 traits). The two added traits pass the heritability gate (z 21.8 and 5.8). Report s for the six clusters in both.
Reading rule (fixed now): a cluster "stays specific" if s >= 2 in the augmented panel. If fewer than half of the six stay specific, the report says the pattern is shared with pain and that the specificity claim is relative to a panel without pain.
No cluster is added or removed from the headline whatever the result.

## P2 Two-machine concordance
Question: do the Mac and Windows pipelines give the same ME/CFS cluster results from the same input? Method: Spearman correlation across the 461 clusters between the Mac model-A p-values (MAGMA, macOS build) and the Windows ones
(`results/clusters_decodeme_gwas_1.tsv`, column p_A), and the same for the specificity score s. Pass rule: Spearman >= 0.95 for -log10 p; a lower value is reported as a reproducibility problem and investigated before anything else is claimed.
