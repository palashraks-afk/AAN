# What this project could discover, and how it will be decided

Status: **NOT YET DETERMINED.** No ME/CFS cell-type result has been opened. The rule (HANDOFF.md, PREREG_v1 G3) is that
the control traits must reproduce known biology first. Everything below was written before the results, so it cannot be
bent to fit them. The final section is filled in only from `results/`.

## 1. Where the project starts: what is already known (so it is not claimed as new)
- ME/CFS has a modest inherited component (liability h2 about 0.095) and 8 genome-wide significant loci (DecodeME, preprint).
- The DecodeME team reports enrichment in brain tissue, and a conference talk and a forum analysis say the signal sits in
  neurons, with medium spiny neurons mentioned as a top cell type. So "the signal is in neurons" is a replication here.
- Independent imaging reports changes in basal ganglia, midbrain, pons, thalamus, hippocampus, amygdala and brainstem.

## 2. The candidate discoveries, each with the rule that decides it
Each is declared only if it passes its pre-registered test, is "confirmed" under the v1.3 rule where that applies
(main test plus at least two of: S-LDSC, permutation null, locus/chromosome drop, independent GWAS), and is not explained
by the generic signal that every brain trait shows (specificity against the 23-trait panel).

| ID | Possible discovery | Decided by | Why it would be new | If it fails |
|---|---|---|---|---|
| A | ME/CFS risk is concentrated in a specific kind of neuron beyond what any polygenic brain trait shows | T2: FDR < 0.05, model B p < 0.05, specificity score s >= 2 against the panel (PREREG v1) | prior work reports "neurons" but not a test against other traits | "The ME/CFS cellular signal is indistinguishable from the generic neuronal signal of polygenic brain traits at this sample size" |
| B | How many independent cell-type axes there are (one axis or several) | L1 forward selection (PREREG v1.5) | **partly covered**: Maccallini 2026 ran FUMA forward selection (one distributed eMSN signal). Ours compares with how many axes schizophrenia, Alzheimer's and RA have | one axis, reported as such |
| C | The signal is brain-restricted rather than body-wide | D2: best neuronal type beats best non-neural type at FDR < 0.05 in gwas_1 and gwas_2 (PREREG v1.2) | **mostly covered at tissue level** (Maccallini: GTEx, no peripheral tissue significant); only the whole-body cell-type resolution is new | "not brain-specific", for example immune or muscle cells also enriched |
| D | The genetic map agrees with independent imaging beyond generic brain traits | E1: literature region set (fixed before results) above the 90th percentile of the panel and permutation p < 0.05 (PREREG v1.3) | convergence between genetics and imaging has not been tested | "no convergence beyond generic brain traits" |
| E | A biological programme sits behind the cell-type signal | D3: Reactome pathways still FDR < 0.05 after conditioning on the top cell type and replicated in gwas_2 (PREREG v1.2) | pathway analysis conditioned on the cell type is not in earlier ME/CFS work | none survive |
| F | Infection-triggered and non-infection ME/CFS differ at the cell-type level | T5 contrast (descriptive; overlapping samples, so no inferential claim) | subtype contrast at cell-type level not found | no difference detectable |
| G | ME/CFS shares its cellular profile with specific other conditions | D4 comorbidity map against the spread of correlations among other trait pairs (PREREG v1.2, v1.6) | cellular rather than gene-level comorbidity | nothing stands out |
| H | Candidate genes ranked by cell-type specificity and genetic support | X1 to X3: only interpreted if the benchmark on migraine, rheumatoid arthritis and IBD passes (PREREG v1) | cell-type-restricted, benchmark-validated ranking | labelled exploratory, no target interpretation |

Full novelty evidence (closest prior work, searches run, grade conditions): `docs/11_novelty_evidence.md`.

## 3. How "unique" will be judged
The uniqueness claim rests on the combination: a pre-registered, panel-calibrated, multi-method confirmed, conditioned
analysis with a benchmark and an independent-literature convergence test, applied to a neglected neurological condition.
No single layer is claimed as new unless the 17 October prior-art search still finds nothing. Each row above will be
written in the report as "to my knowledge not previously reported (searched [date])", never as "first".

## 4. Final statement (filled from the results, not before)
- Gate G3 (controls): **PASS** (2026-10-09 19:36). Schizophrenia: excitatory-neuron cluster at FDR < 0.05. Alzheimer's: microglia (FDR 3e-7). Height: strongest three clusters all fibroblasts, none neuronal. Output in `results/G3_check.txt`.
- ME/CFS localisation (T1), **preliminary, first look 2026-10-09, before calibration and confirmation**: 11 of 461 clusters at FDR < 0.05 in model A (MAGMA gene-property, DecodeME gwas_1): deep-layer intratelencephalic (6), amygdala excitatory (2), upper-layer intratelencephalic (2) and one splatter cluster; strongest z about 3.6 (smallest p 1.7e-4). In model B (also conditioning on mean neuronal expression) no cluster reaches FDR < 0.05 (best FDR 0.067). The medium spiny neuron superclusters are not significant here (best z 2.8). The signal is modest compared with the controls (schizophrenia z about 7, Alzheimer's about 6). Reading: a weak glutamatergic (excitatory) cortical and amygdala signal that has not yet been tested for specificity. Not a discovery until T2 (panel calibration), confirmation (C1 to C4) and gwas_2 replication are run.
- Discovery A to H: [pending, one line each: declared / not declared, with the number behind it and the checks passed]
- The single most defensible new statement: [pending]
- What the results do NOT show: [pending, always includes: association not causation; three-donor atlas; CNS only; European
  ancestry; hypotheses for researchers, not treatment advice]
