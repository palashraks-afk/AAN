# Extension layers (Mac pipeline) — running alongside the main pipeline

Started 2026-10-10. This file is separate from `HANDOFF.md` and the Windows run: it records extra analysis layers run on the Mac
after the control gate G3 passed on the Windows machine. The main pipeline (`run_everything.py`) is not changed by it.
Every layer is pre-registered (`preregistration/PREREG_v1.7_specific_component.md`, tag `prereg-v1.7`) before its ME/CFS result was run.
**Status of results: see the table; nothing is called a discovery until its pre-registered rule and checks are met.**

## Why these layers
The main pipeline asks which cell types carry ME/CFS signal and whether that beats other traits cluster by cluster. Prior work
(Maccallini et al. 2026) already reports cell-type enrichment, so localisation is a replication. These layers ask what is left over:

| Layer | Question | Method | Decided by | Status |
|---|---|---|---|---|
| V (gate) | Does the new method work on known biology? | OLS plus permutation vs MAGMA on controls | PREREG v1.7 V1 to V3 | **failed** (schizophrenia agreement 0.768 < 0.8; height check stricter than G3). Recorded in `results/x_R1_validation.tsv`; no ME/CFS output of this version was run |
| R1' (v1.9) | Same question, done inside MAGMA: condition on the principal components of the other panel traits | MAGMA gene-property conditioning on avg_all, avg_neuron, PC1 to PC5 (leave-one-out panel for the s score) | FDR < 0.05 gwas_1, p < 0.05 gwas_2, s >= 2 | **done: not declared** (gate V' passed; best FDR 0.073) |
| R1 (v1.7, superseded) | What cell types carry the ME/CFS signal that is **not** shared with 20+ other traits? | gene-level z of ME/CFS residualised on the principal components of the panel's gene z, then cell-type regression with a size-and-expression-matched permutation null | FDR < 0.05 in gwas_1, p < 0.05 in gwas_2, s >= 2 | pending |
| R2 | Which genes carry that specific component? | residual z >= 3.5, overlap with DecodeME loci and approved-drug targets | descriptive; benchmark rule before any target is read as a hypothesis | pending |
| F1 | Do drug-relevant gene families carry signal beyond brain expression? | MAGMA competitive set test on 8 pre-fixed HGNC families, conditioned on expression | FDR < 0.05 gwas_1 and p < 0.05 gwas_2 | **done: not supported.** Best are ion channels (p 0.021, FDR 0.12) and nuclear receptors (p 0.030, FDR 0.12) in gwas_1, same direction and similar size in gwas_2 (p 0.027, 0.041; FDR 0.16). None passes FDR < 0.05. |

## Files
- `scripts/x_specific_component.py`: V, R1, R2 (`--validate`, `--mecfs`)
- `scripts/sync_mac.sh`: builds `~/AAN_mac` (macOS paths) from this repo so commits keep the Windows paths
- outputs: `results/x_*.tsv` (R1 profiles, calls, validation gate, top residual genes)

## Rules this extension follows
1. Validation gate V first; if it fails, R1 is reported as a failed method and not interpreted.
2. Calls use the same style of rule as the main pre-registration: FDR, replication in gwas_2, and a specificity score against the panel.
3. A null ("nothing specific remains beyond what the panel traits share") is a reportable result and is pre-committed to.
4. No "first", "never done", "cure" or "treatment" wording. Novelty is graded in `SCORES.md`, not here.
5. Mac and Windows results are not mixed in one table; the Mac panel results are used only by these layers.

## Results
### F1 (2026-10-10), gene families: observed, not supported
Eight pre-fixed HGNC families, MAGMA competitive test conditioned on mean expression and mean neuronal expression, BH across 8.
Ion channels (183 genes) and nuclear receptors (21 genes) show nominal signal in both ME/CFS GWAS (gwas_1 p = 0.021 and 0.030; gwas_2
p = 0.027 and 0.041) but FDR is 0.12 to 0.16, so by the pre-registered rule neither is declared. gwas_2 shares all cases with gwas_1, so
it is a robustness check, not an independent replication. As a reference, Alzheimer's shows no channel signal; its best family is
immunoglobulin-related genes (p 0.015, FDR 0.12). Output: `results/x_F1_families.tsv`. Reading: a hint worth following up with a larger
GWAS, not a finding.

### R1' (2026-10-10), signal left after removing what 19 other traits share: not declared
Panel of 19 traits (insomnia dropped by the heritability gate, h2 z = 1.4). Gate V' passed first: Alzheimer's keeps a microglia cluster at
FDR < 0.05 under leave-one-out conditioning, and height's three strongest clusters are all fibroblasts. For ME/CFS gwas_1, conditioning on
the first five principal components of the panel's gene-level z leaves **no cluster at FDR < 0.05** (best 0.073), so by the pre-registered
rule nothing is called ME/CFS-specific by this layer. The clusters at the top are the same ones the main pipeline found: upper- and
deep-layer intratelencephalic neurons (ULIT_121 z 3.19, DLIT_150/151/152 z 3.1 to 3.2), the amygdala excitatory cluster Amex_153 (z 2.95,
s 3.69, the highest specificity score), and two eccentric medium spiny neuron clusters (EMSN_224 z 2.95, s 3.51; EMSN_222 z 2.79, s 2.78).
All pass p < 0.05 in gwas_2 (same cases, different controls). Reading: the amygdala and striatal-type clusters keep the highest
specificity scores, but their signal is too weak to clear an FDR bar once the shared genetics is conditioned out. This agrees with the
main pipeline's own caveat (none pass FDR under model B) and does not strengthen the amygdala claim to "confirmed". It is
consistent with, and does not add to, the earlier reports of medium spiny neurons. Outputs: `results/x_R1p_mecfs.tsv`, `results/x_R1p_validation.tsv`.

### GTEx check (v1.10, 2026-10-10): the amygdala result does not replicate in independent bulk data
Pre-registered H1 (Brain - Amygdala vs other brain regions, within-brain specificity): p = 0.49 (gwas_2 p = 0.35), specificity 0.25: **not supported**.
H2 (Nerve - Tibial, a peripheral nerve): p = 0.75, not supported. H3 (Pituitary): p = 7.6e-4 and gwas_2 p = 0.018, but specificity against 19 traits is
1.1 (< 2), so **not supported** (pituitary is enriched, but not more than other traits). Across all 54 tissues every brain tissue is enriched
(FDR about 0.001, caudate and nucleus accumbens strongest) with specificity below 1, which only confirms DecodeME's brain enrichment. Bulk tissue
cannot resolve a small neuronal population, so this weakens but does not refute the single-nucleus amygdala result. Outputs: `results/x_GTEx_S54.tsv`, `results/x_GTEx_SB.tsv`.
Deep novelty check of the whole claim: `docs/12_discovery_novelty_check.md`.

### R2: not run (it depended on the OLS version of R1).

### Independent single-nucleus atlas (PREREG v1.11 P1, Tran et al. 2021 amygdala + nucleus accumbens), 2026-10-10: amygdala excitatory neurons not replicated (weak test)
35 cell-type groups (15 amygdala, 20 nucleus accumbens; 5 and 8 donors), MAGMA conditioned on mean expression, specificity against the 19-trait panel.
- **H-A1 (pre-registered): not supported.** AmyG Excit_A p = 0.12 (gwas_2 p = 0.079, s = 0.90) and AmyG Excit_C p = 0.15 (gwas_2 p = 0.13, s = 1.48); the bar was p < 0.025 and s >= 2.
- Nothing reaches FDR < 0.05 (best FDR 0.20). The nominal leaders are inhibitory amygdala groups (Inhib_H p = 0.006, gwas_2 p = 0.024; Inhib_A p = 0.017, s = 1.97) and
  indirect-pathway medium spiny neurons of the nucleus accumbens (MSN.D2_B p = 0.013, gwas_2 p = 0.006; D2_D and D1_D/E at p of 0.04 to 0.08 in gwas_2). Exploratory, not supported.
- **Why this is a weak test:** the dataset has only 443 glutamatergic amygdala nuclei (Excit_C has 55, Excit_B 44 and is excluded), the labels are coarse, and the Siletti
  clusters that carried the signal (Amex_153/175) have no one-to-one match here. A null is uninformative about the Siletti result. Output: `results/x_ATLAS_tran.tsv`.
- The direction of the MSN.D2 signal is consistent with earlier informal and Maccallini reports of medium spiny neurons, and is the only nominal pattern that also reaches p < 0.05 in gwas_2, but it is not corrected-significant.
