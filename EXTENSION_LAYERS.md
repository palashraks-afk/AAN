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
| V (gate) | Does the new method work on known biology? | OLS plus permutation vs MAGMA on controls; Alzheimer's keeps microglia after residualisation; height stays non-neuronal | PREREG v1.7 V1 to V3 | pending (waiting for panel MAGMA) |
| R1 | What cell types carry the ME/CFS signal that is **not** shared with 20+ other traits? | gene-level z of ME/CFS residualised on the principal components of the panel's gene z, then cell-type regression with a size-and-expression-matched permutation null | FDR < 0.05 in gwas_1, p < 0.05 in gwas_2, s >= 2 | pending |
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

### V, R1, R2: waiting for the comparison panel's MAGMA runs (the panel needs about 20 traits; downloads and MAGMA are still running).
