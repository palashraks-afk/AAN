# What this project could discover, and how it will be decided

Status: **results in (2026-10-10), see section 4.** The control traits reproduced known biology first (G3 PASS). Sections 1 to 3
were written before the results; section 4 is filled only from `results/`.

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
- Full pipeline finished 2026-10-10 (`results/RUN_REPORT.md`, `results/EVIDENCE_LEDGER.md`, `results/layers_summary.txt`).
  All 11 clusters from T1 are "confirmed" under the v1.3 rule (S-LDSC, permutation null with empirical p <= 0.002, locus and
  chromosome drop, gwas_2). **Caveat for C4:** gwas_2 uses the same cases with a different control subset (cluster z correlation
  0.978 with gwas_1), so it tests stability to the control set and is not an independent replication.
- Discovery A to H, decided by the rules in section 2:
  - **A, declared, modest.** Six of the 11 clusters pass FDR < 0.05, model B p < 0.05 and s >= 2 against the panel:
    amygdala excitatory Amex_153 (s 2.98) and Amex_175 (2.50); deep-layer intratelencephalic DLIT_152 (2.03) and DLIT_150 (2.01);
    upper-layer intratelencephalic ULIT_121 (2.35); splatter Splat_402 (2.41). The other five (DLIT_136, 147, 148, 151, ULIT_126)
    are real but shared with other brain traits. Limits: z is 3.0 to 3.4 and FDR 0.03 to 0.05; with the stricter model B (FDR) none pass
    (best 0.067); with model C (also conditioning on depression, BMI and insomnia genetics) 7 clusters pass FDR < 0.05, of which 4
    are among the six (DLIT_152, DLIT_150, ULIT_121, Amex_153); cluster-level lambda is 1.48; the panel has 19 traits, so s is
    noisy. Amygdala excitatory neurons are the most specific.
  - **B, declared as "one axis".** Forward selection keeps one cluster (DLIT_136) for ME/CFS, Alzheimer's and RA, two for
    schizophrenia. Same conclusion as Maccallini 2026, so a replication.
  - **C, not declared.** 12 of 154 whole-body cell types are at FDR < 0.05, including non-neuronal ones (somatotrophs, lactotrophs,
    oligodendrocyte precursors, oligodendrocytes). Only somatotrophs (s 2.08) and one other are "supported" by gwas_2.
  - **D, not declared.** Literature-region convergence statistic -0.455, 10th percentile of the panel, permutation p 0.91.
  - **E, not declared.** No Reactome pathway reaches FDR < 0.05 (best 0.69).
  - **F, descriptive only (overlapping samples).** Clusters at FDR < 0.05: infectious onset 55 (12 specific, top DLIT_146),
    non-infectious onset 11 (11 specific, top MGE_259, an interneuron cluster), female 0, male 0. The counts mostly reflect sample
    size and power; no claim of a difference.
  - **G, descriptive.** Highest profile correlations with ME/CFS: education years (0.78), schizophrenia (0.77), BMI (0.60),
    neuroticism (0.59). They sit above the 90th percentile of other trait pairs, but this mostly reflects the shared neuronal profile.
  - **H, exploratory.** The benchmark passed for 1 of 3 diseases (IBD; migraine and RA failed), and the rule needs 3, so gene
    rankings (top: ISL1, CACNA1E, STT3B, DCC, PCDH17) are labelled exploratory with no target interpretation.
- **Robustness of A to the choice of panel (PREREG v1.8 R2 to R4, run 2026-10-10, `results/panel_robustness.tsv`).** Rule fixed beforehand: robust
  if s >= 2 in at least 90% of leave-one-trait-out panels and at least 80% of random half panels. **Only the two amygdala excitatory
  clusters pass**: Amex_153 (s 2.89 to 4.02 across every variant) and Amex_175 (s >= 2.4; 88% of half panels). The other four are
  panel-dependent: ULIT_121 (77% of half panels), Splat_402 (80%, borderline), DLIT_150 (37% of leave-one-out, s about 1.94 to 2.0) and DLIT_152
  (53%, s about 1.96 to 2.0). Updated reading of A: the specific, panel-robust result is **amygdala excitatory neurons**; the cortical
  clusters are borderline and should be described as such.
- **R1 (v1.8) permutation null at the pre-registered 5,000 permutations (2026-10-10):** all 11 clusters keep empirical p <= 0.0012
  (nine at the minimum 0.0002, Splat_402 0.0012). The 1,000 run is kept as `results/permutation_decodeme_gwas_1_n1000.tsv`. This
  confirms the C2 check; it says nothing new about specificity.
- **Method validation (PREREG v1.12, run 2026-10-10, `results/panel_as_target.tsv`).** Each panel trait was made the target with the other traits as its
  panel. Rule fixed beforehand: the score recovers the expected class in at least 4 of 5 (Alzheimer's microglia, RA and IBD immune/vascular, schizophrenia
  excitatory neurons, height non-neuronal). **Result: 4 of 5 recovered** (Alzheimer's, RA, IBD, height yes; schizophrenia no: its neuronal signal is shared with
  other brain traits, so none of its 280 clusters is specific). Among traits with at least 3 significant clusters (9 traits), the median share of significant
  clusters that are specific is 3%; the heavily neuronal traits are at 0% (schizophrenia), 0.5% (education), 1.3% (neuroticism) and 3% (BMI). ME/CFS is 6 of 11 (55%),
  above all four of those but **below the 90th percentile (70%)** of traits, so by the rule nothing "unusual" is claimed. Caveat: ME/CFS has few significant clusters, and
  a small count makes the share noisy and tends to favour the strongest tail. What it supports: the score separates traits whose signal is generic (neuronal brain traits)
  from traits with a distinct class (immune, connective tissue), and ME/CFS looks more like the second group than the first.
- **Solidity checks W1 to W8 (PREREG v1.14, run 2026-10-10, `results/w_checks.txt`; W2, W9, W10 are still running).**
  - W1 matched-power panel: all six keep s >= 2 when the panel is only the 13 low-power traits (Amex_153 s 4.01, Amex_175 2.59, DLIT_152 2.44, DLIT_150 2.32, ULIT_121 2.18,
    Splat_402 2.12). Passes. (CRP has 34 significant clusters, so by the stated cutoff of 30 it counts as a strong trait.)
  - W3 chromosome jackknife of s: the intervals are wide and **no cluster meets the "stable" rule (lower bound > 1.5)**: Amex_153 0.58 to 5.21, Amex_175 1.07 to 3.79, the others
    0.4 to 4.1. The point values hold (every leave-one-chromosome-out s is at least 1.65, Amex_153 at least 2.51). Reading: the specificity score is imprecise at this sample size.
  - W4 **where the cells come from: Amex_153 is 93% amygdala dissections, but Amex_175 is 91% cerebral cortex dissections (1% amygdala).** The "amygdala excitatory" name of Amex_175 is a
    supercluster label; its cells are mostly cortical. **The clean amygdala claim therefore rests on Amex_153 alone.** Reword: "an excitatory cluster from amygdala dissections (Amex_153) and a
    related cortical cluster".
  - W5 subsets: all six clusters pass the Bonferroni threshold (p < 0.0021) in the infection-onset subset (e.g. Amex_153 p 6.8e-5); none do in non-infectious onset, female-only or male-only.
    Descriptive (overlapping samples, and the subsets differ in size); it suggests the signal sits mainly in infection-triggered ME/CFS.
  - W6 amygdala region score: rank 2 of 17 regions for ME/CFS (Paleocortex first) but the 79th percentile of the panel, **below the pre-set 90th, so the rule fails**.
  - W7 power projection (an assumption, not a result): with the 21,620 recruited cases (z scaled by 1.18), model B would have 24 clusters at FDR < 0.05 and all six would pass (best FDR 0.009 to 0.021),
    and model C 42. If the effect is real, the final DecodeME data should settle it.
  - W8 other definitions: all six are definition-robust (rank among the panel at the top for all 19 traits, weighted s 2.17 to 3.57, and top 10% for ME/CFS but for at most 4 of 19 panel traits).
  - W9 drivers (`results/w9_drivers.tsv`): dropping the top 1, 5 and 20 contributing genes lowers z gradually (Amex_153 3.43, 3.40, 3.25, 2.81; Amex_175 3.16, 3.10, 2.91, 2.36) and both clusters stay
    nominally significant after dropping 20 genes (p 0.0025 and 0.0092). **Passes: the signal is spread over many genes, not a few.** Top contributors for Amex_153 include 57554 (LRRC7).
- **Later checks that did not strengthen A (run by the second session, 2026-10-10; see `EXTENSION_LAYERS.md`, `docs/12_discovery_novelty_check.md`).**
  (1) R1' (PREREG v1.9): conditioning inside MAGMA on the first five principal components of the 19 other traits leaves **no cluster at
  FDR < 0.05** (best 0.073); Amex_153 still has the highest specificity score (3.69). (2) GTEx independent check (v1.10, bulk tissue,
  different donors): Brain - Amygdala is significant against all tissues like every brain tissue, but is **not special among brain regions**
  (p 0.49). Bulk tissue may hide a small cell population, so this does not refute the single-cell result, but it gives no support either.
  Reading: A meets its pre-registered rule, but the evidence is modest and fragile. It is to be written as a hypothesis, never as a
  confirmed discovery.
- **Prior-art check 2026-10-10 (web search; to be repeated in PubMed, bioRxiv, medRxiv on 17 Oct).** Not new as localisation: an informal
  forum analysis of the same Siletti atlas (S4ME thread "Using MAGMA on ME/CFS genetic data": forestglip, Tralfamadorian97, ME/CFS
  Science Blog) already noted amygdala excitatory neurons, deep-layer and upper-layer intratelencephalic cells next to eccentric medium
  spiny neurons, and Maccallini's blog (4 Oct 2025) lists amygdala among enriched regions and layer 2/3, 5, 6 glutamatergic neurons.
  Those analyses are uncorrected for other traits and are not peer reviewed. Not found anywhere: a test of whether these clusters are
  specific to ME/CFS against a panel of other traits, with conditioning, permutation null and locus drop. So the new part is the
  calibration, not the cell types.
- The single most defensible new statement: earlier informal analyses pointed at amygdala and intratelencephalic excitatory neurons;
  here a calibrated, pre-registered test against 19 other traits shows that a small subset of them (mainly amygdala excitatory
  neurons) carries ME/CFS signal beyond what generic brain traits show, while the broader neuronal enrichment is shared. A modest
  result: the cell types are not new, the specificity test is.
- What the results do NOT show: association not causation; three-donor atlas; CNS only; European ancestry; hypotheses for
  researchers, not treatment advice; gwas_2 is not an independent replication; the evidence is modest (FDR 0.03 to 0.05, none
  under model B FDR); the medium spiny neuron finding of earlier work is not reproduced here.
