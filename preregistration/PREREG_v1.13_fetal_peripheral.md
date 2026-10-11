# Pre-registration v1.13: peripheral and autonomic nervous system cell types, and fetal limbic neurons (Cao et al. 2020 atlas)

Written 2026-10-10 (Mac session, task M-extra), tag `prereg-v1.13`. Adds analyses only. No ME/CFS result on this atlas has been computed; I inspected only the cell-type labels and counts.

## Why
Every brain atlas used so far lacks peripheral nerves and autonomic neurons, which the project's own limitations section names and which matter in ME/CFS (autonomic
dysfunction; Lee et al. 2026 saw enteric nervous system neurons in a fetal atlas without any control for other traits). The Cao et al. 2020 human fetal atlas
(CELLxGENE collection c114c20f-1ef4-49a5-9c2e-d965787fb90c, the 1-million-cell subset, 10 to 18 weeks of gestation, 15 organs, 77 `Main_cluster_name` types) contains
peripheral neuron types and fetal limbic neurons. It also gives a developmental-timing angle.

## Method
Raw counts aggregated per `Main_cluster_name` (groups with at least 50 nuclei, all kept), specificity as in v1.10/v1.11 (fraction of expression across groups, rank-normalised),
MAGMA gene-property conditioning on mean expression, one-sided positive, gwas_1 and gwas_2, specificity score s against the same 19 panel traits.

## Hypotheses fixed now
- H-P1 (peripheral nervous system): six groups: Sympathoblasts, ENS neurons, Visceral neurons, Schwann cells, ENS glia, Chromaffin cells. "Supported" if, for any of them, p < 0.05/6 in gwas_1,
  p < 0.05 in gwas_2 (same direction) and s >= 2.
- H-P2 (fetal limbic neurons): the single group "Limbic system neurons". "Supported" if p < 0.05 in gwas_1 and in gwas_2 and s >= 2.
- Exploratory: all 77 groups with BH-FDR, no confirmation language.
- Stated limitation: this atlas is a 1-million-cell subsample of 4 million, the peripheral groups have only a few hundred nuclei, labels are coarse, and the cells are fetal. A null is weak evidence.

## Not allowed
Adding or dropping groups or changing thresholds after any ME/CFS result on this atlas is seen.
