# Pre-registration v1.17: amygdala excitatory neurons against cortical and hippocampal excitatory neurons (task M2 amendment to v1.11 P1)

Written 2026-10-10 (Mac session), tag `prereg-v1.17`, before the three extra regions were opened or any ME/CFS result on them computed.
**This is a post-hoc amendment.** v1.11 P1 used only the amygdala and nucleus accumbens datasets of Tran et al. 2021 and gave a null for the amygdala excitatory groups (p = 0.12, 0.15).
It left the question open because that atlas had no cortical or hippocampal excitatory neurons to compare against, so "amygdala excitatory" could not be told apart from "excitatory". I am adding
the other three regions of the same atlas (dorsolateral prefrontal cortex, subgenual anterior cingulate cortex, hippocampus). The v1.11 result stays reported as it was. This extension is exploratory and cannot be called a confirmation.

## Method
Same as v1.11 P1: raw counts aggregated per `author_cell_type` per region, groups with >= 50 nuclei, `drop.*` excluded, specificity across all groups of all five regions, MAGMA conditioning on mean expression,
one-sided positive, gwas_1 and gwas_2, specificity score s against the same 19 traits.

## Hypothesis fixed now
H-A3: AmyG Excit_A (344 nuclei) or AmyG Excit_C (55 nuclei) has p < 0.025 in gwas_1 and p < 0.05 in gwas_2 with s >= 2 in the five-region analysis. All other groups are reported with BH-FDR as exploratory.
Also reported: the best cortical/hippocampal excitatory groups, so the reader can see whether the signal belongs to excitatory neurons in general.
Stated: a null remains weak evidence because of the small number of amygdala excitatory nuclei.

## Not allowed
Adding regions or groups, or changing thresholds, after any result of this extension is seen.
