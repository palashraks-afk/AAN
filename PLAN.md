# PLAN — ME/CFS nervous-system genetic atlas + treatment-hypothesis dossier

Written 2026-10-07. Deadline 2026-10-20 23:59 CT. Status: plan only, nothing run.

---------------------------------------------------------------------------------------------------
## 0. One-paragraph summary

ME/CFS (myalgic encephalomyelitis / chronic fatigue syndrome) affects over a million people in the US
alone, is classified by WHO as a neurological disease (ICD-10 G93.3), and is among the most
under-studied large conditions. DecodeME (University of Edinburgh) released the largest ME/CFS GWAS
(15,579 cases, 84% female, European ancestry; 259,909 UK Biobank controls) with open summary
statistics. The team reported brain-tissue enrichment, but I found no calibrated nervous-system
**cell-type** analysis. This project (1) maps the signal onto human brain / brainstem / hypothalamus
cell types (Siletti 2023 Human Brain Cell Atlas), (2) tests whether any cell type is enriched
**beyond what any polygenic brain trait shows** (specificity calibration against a panel of other
traits), (3) tests whether the signal survives removing shared genetics with depression, BMI and
insomnia, (4) compares infection-triggered vs non-infection-triggered and female vs all cases, and
(5) turns genes driving the specific signal into a direction-aware, safety-checked,
benchmark-validated **drug-target dossier** (hypotheses only). Everything is pre-registered.

---------------------------------------------------------------------------------------------------
## 1. Why this idea (and what it is NOT)

- Large: ME/CFS is far above the 50,000-people bar (>1 million US; ~17 million worldwide is often cited
  — verify the source before quoting).
- Neglected: few studies, stigma, no approved disease-modifying treatment.
- Human genetics: open summary statistics, no human-subjects data handling by the applicant.
- Neuroscience relevance: the question is *which cells of the nervous system* carry the genetic signal.
- NOT claimed: a cure, a clinical recommendation, or "first ever". Honest grades are in
  docs/04_honest_grades.md. Novelty is graded 4, not 5.

---------------------------------------------------------------------------------------------------
## 2. Research question and hypotheses

Primary question: In which human nervous-system cell types does the genetic risk for ME/CFS act, and
is that pattern specific to ME/CFS?

- **H1 (localisation).** ME/CFS GWAS signal is enriched in specific brain / brainstem / hypothalamic
  cell types at FDR < 0.05 across clusters (thresholds in PREREG).
- **H2 (specificity).** Those enrichments exceed the null distribution formed by a panel of ≥ 20
  other traits (i.e. they are not the generic "neurons" signal of any polygenic brain trait).
- **H3 (not just mood/body-size).** Enrichment persists after conditioning out shared genetic signal
  with depression, BMI and insomnia.
- **H4 (subtype).** The infection-triggered subset differs in cell-type profile from the
  non-infection subset; female-only resembles all-cases. (Low power; descriptive unless gated in.)
- **H5 (translational).** Genes driving the specific cell-type signal can be ranked into evidence
  tiers (colocalisation + causal direction + existing drug with matching direction + safety), and the
  ranking pipeline recovers known successful drug targets in positive-control diseases better than
  chance.

A null result (no ME/CFS-specific cell type) is a valid, reportable finding and is pre-committed to.

---------------------------------------------------------------------------------------------------
## 3. Data (all open; see data/README.md for sizes and MD5)

### 3.1 ME/CFS GWAS — DecodeME, OSF project rgqs3, folder "DecodeME Summary Statistics"
REGENIE output, GRCh38, 8,835,520 variants, European ancestry, controls from UK Biobank.

| Analysis | Cases | Controls | Use |
|---|---|---|---|
| gwas_1 (primary) | 15,579 | 259,909 | main enrichment |
| gwas_2 | 15,579 | 155,790 (a different UKB control subset) | internal replication of enrichment |
| gwas_1_female | 12,833 | 218,949 | sex check |
| gwas_1_male | 2,746 | 40,960 | probably too small — gated out by h2 z-score |
| gwas_1_infectious_onset | n/a (subset) | UKB | subtype contrast (H4) |
| gwas_1_non_infectious_onset | n/a (subset) | UKB | subtype contrast (H4) |
| gwas_qced.var.gz, imputed.info.gz | — | — | variant QC / INFO scores |

Case definition: Canadian Consensus Criteria-based inclusion/exclusion. SNP-heritability (liability
scale) reported by the team: h2 = 0.095 (SD 0.006) from GWAS-1.
NOTE: GWAS-1 and GWAS-2 share all cases; female/male and infection/non-infection subsets overlap
GWAS-1 — they are NOT independent replications. State this in the report.

### 3.2 Cell-type reference — Human Brain Cell Atlas v1.0 (Siletti et al. 2023, Science)
~3 million nuclei, ~100 dissections across forebrain, midbrain, hindbrain; 31 superclusters,
461 clusters, 3,313 subclusters. Use the **aggregated** genes × clusters loom
(`adult_human_20221007_Pool_Clean.agg.loom`), NOT the 37 GB full loom.
Source: github.com/linnarsson-lab/adult-human-brain (BSD-2 code; data licence not separately stated).
Check file size with an HTTP HEAD request before downloading.

### 3.3 Comparison-trait panel (specificity null) — availability UNVERIFIED, check each
Aim ≥ 20 traits, mixing brain-relevant and non-brain controls. Candidates:
schizophrenia (PGC3) [positive control: cortical neurons], Alzheimer's disease (Bellenguez 2022)
[positive control: microglia], height [negative control: non-neural], depression (PGC),
insomnia, BMI/obesity (GIANT), neuroticism, migraine, multiple sclerosis (IMSGC), chronic pain,
fibromyalgia (2026), long COVID (COVID-19 HGI / 23andMe), Parkinson's, ADHD, bipolar, educational
attainment, CRP, IBS, autoimmune controls (RA, T1D) [immune-cell positive control].
Record for each: source URL, N, build, licence, and whether UK Biobank controls overlap with DecodeME
controls (matters for any cross-trait genetic-correlation/conditioning step).

### 3.4 Translational resources (open) — UNVERIFIED access, check each
Brain / single-nucleus eQTL (GTEx v10 brain; SingleBrain snRNA-eQTL meta-analysis), pQTL (UKB-PPP;
brain pQTL from ADSP FunGen xQTL — large), Open Targets Platform (GraphQL API), DGIdb, ChEMBL,
gnomAD constraint, FinnGen R13 LoF burden browser.

---------------------------------------------------------------------------------------------------
## 4. Methods

### 4.1 Harmonise and QC (day 2)
- Read REGENIE files; columns expected: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE
  CHISQ LOG10P EXTRA (VERIFY against the README_shared_sum_stats.txt — I could not read it; OSF
  returned HTTP 429).
- Convert LOG10P → P. Filter INFO ≥ 0.8, MAF ≥ 0.01, drop MHC (chr6:25–34 Mb, GRCh37 — use the
  GRCh38 equivalent) for LDSC.
- Genome build: DecodeME is GRCh38; most LD reference panels (1000G Phase 3) are GRCh37. Decide ONE
  route and record it: lift variants over (pyliftover), or use a GRCh38 LD reference. Do not mix.

### 4.2 Heritability gate (day 2)
LDSC SNP-heritability and z-score for every GWAS. A GWAS is analysed only if h2 z ≥ 4 (PREREG).
Expect gwas_1_male to fail. This is a pre-committed exclusion, not p-hacking.

### 4.3 Cell-type specificity matrix (day 3)
From the aggregated loom compute, per cluster, a gene specificity score (e.g. expression-weighted
mean / fraction of total across clusters, CELLEX-style ES_mu). Per cluster take the top 10% most
specific genes (and run 5% / 20% as sensitivity). Control for general neuronal expression by
including mean expression across all clusters as a covariate in the gene-property model.

### 4.4 Enrichment (days 4–5)
- Gene-level association: MAGMA v1.10+ (older versions have a documented type-I error issue in
  SNP→gene aggregation). Window 35 kb up / 10 kb down; 1000G EUR reference.
- Gene-property (continuous specificity) and competitive gene-set tests per cluster.
- Orthogonal confirmation with stratified LDSC (LDSC-SEG) using the same gene sets.
- Multiple testing: Benjamini–Hochberg across clusters within a trait; report both cluster and
  supercluster level.
- Run the **identical** pipeline on every panel trait.

### 4.5 Specificity calibration (day 5) — the main novelty
For each cluster, compare the ME/CFS enrichment statistic with the empirical distribution of that
cluster's statistic across the panel traits. Report (i) raw enrichment, (ii) rank/percentile among
traits, (iii) enrichment conditional on the generic-neuron factor. A cluster counts as
"ME/CFS-specific" only if it passes all three rules frozen in PREREG.

### 4.6 Conditioning on shared genetics (day 6)
Preferred: gene-level conditional analysis — re-run MAGMA gene-set tests with the gene-level
association of depression / BMI / insomnia as covariates (MAGMA conditional model). Fallback: drop
genes inside loci significant for those traits and re-test. mtCOJO/CTG-COJO only if sample overlap
is handled (UK Biobank controls are shared between DecodeME and many UKB-derived comparison GWAS —
this inflates cross-trait correlations; document and prefer the gene-level approach).

### 4.7 Subtype and sex contrasts (day 6)
Female vs all; infection vs non-infection onset. Report effect sizes with CIs. Because the subsets
overlap, use differences in enrichment effect size with a permutation/bootstrap CI, and call it
descriptive.

### 4.8 Translational dossier (days 7–9) — hypotheses only
1. Gene prioritisation: genes ranked by (MAGMA gene Z) × (specificity in the ME/CFS-specific
   clusters). Not "all GWAS hits".
2. Causal evidence: colocalisation (coloc / SuSiE-coloc) of ME/CFS loci with brain eQTL/pQTL.
   R is not installed — either install R + coloc/TwoSampleMR or implement approximate-Bayes-factor
   colocalisation in Python (document which).
3. Direction: Mendelian randomisation on cis-eQTL/pQTL instruments → does *more* gene product raise
   or lower risk → inhibitor vs activator is the therapeutically relevant direction.
4. Druggability: Open Targets, DGIdb, ChEMBL — keep only drug–target pairs whose mechanism direction
   matches step 3.
5. Safety: gnomAD LoF constraint (pLI / LOEUF), LoF carriers in FinnGen/UKB if accessible, breadth
   of expression.
6. Evidence tiers (frozen in PREREG): A = coloc + consistent MR direction + direction-matched
   existing drug + no major constraint flag; B = three of four; C = candidate only.
7. **Benchmark (non-negotiable):** run the identical pipeline on positive-control diseases with
   known successful drug targets (e.g. migraine → CGRP pathway; rheumatoid arthritis → JAK/IL6R;
   multiple sclerosis → known targets). Report how highly known targets rank (AUC / rank
   percentile) versus random gene sets matched for size and length. If the pipeline cannot
   rediscover known winners, the ME/CFS dossier is reported as exploratory and not interpreted.
   Context: genetic support is associated with ~2.6× higher probability of drug success
   (Minikel et al., Nature 2024) — observational, not a guarantee for any target.

### 4.9 Stretch (only if everything above is finished by 14 Oct)
Repeat 4.2–4.5 on FinnGen R13 neglected-condition endpoints (trigeminal neuralgia ~2,356 cases;
cluster headache ~1,600–2,050 cases) and fibromyalgia/long COVID if sumstats are open; apply the
heritability gate; attach each condition's "neglect ratio" (burden vs publications/funding/trial
sites — see docs/01 for the Neglect Atlas design). Drop without regret if time is short.

---------------------------------------------------------------------------------------------------
## 5. Decision rules and fallbacks

| Situation | Action |
|---|---|
| Tools won't install by end of 9 Oct | Use FUMA web server for MAGMA gene/gene-set; do specificity calibration in Python on exported output |
| h2 z < 4 for a subset | Exclude subset; report exclusion |
| No ME/CFS-specific cluster survives calibration | Report as null; the finding is "ME/CFS looks like a generic polygenic neuronal signal" |
| Benchmark fails | Dossier labelled exploratory; no target interpretation |
| Core enrichment not working by 13 Oct | Drop translational/stretch; submit rigorous tissue+cell-type + calibration + limits |
| Someone publishes the same specificity analysis | Cite them, reposition as independent replication, lower novelty claim |

---------------------------------------------------------------------------------------------------
## 6. Deliverables

- Abstract ≤ 300 words (applicant-written; first-round judging is on this alone).
- Research report (figures/tables labelled and readable): intro, methods, results, interpretation,
  pitfalls & limitations, future work.
- Bibliography.
- Public reproducible pipeline in this repo + ranked hypothesis list.
- Optional: send the dossier to the DecodeME team and a patient organisation (e.g. Action for ME)
  and to a clinician for review — the realistic route to real-world impact.

---------------------------------------------------------------------------------------------------
## 7. Known limitations to state openly

- GWAS cases are 84% female; male analysis is underpowered.
- European ancestry only.
- Controls are UK Biobank (healthy-volunteer bias; overlap with comparison GWAS).
- Case definition differs across datasets; DecodeME cases are volunteer-recruited via a questionnaire.
- Cell-type enrichment is correlational; atlas is from 3 postmortem donors.
- Specificity calibration depends on the chosen trait panel; report sensitivity to panel composition.
- Brain atlas covers CNS; peripheral nervous system / immune cells matter in ME/CFS and may be
  under-represented. Say so.
- Genetic support for a target does not mean a drug will work or is safe. No clinical advice.

---------------------------------------------------------------------------------------------------
## 8. Open verifications (do these first)

1. Read `README_shared_sum_stats.txt` (column definitions) — OSF returned 429; retry later.
2. Confirm the three signatures (parent/guardian, teacher, mentor) have been requested.
3. Re-search for prior art on: ME/CFS cell-type enrichment with specificity calibration; ME/CFS
   target prioritisation with direction-aware MR; whether DecodeME's final paper added cell-type work.
4. Confirm availability/licence of every comparison-trait GWAS.
5. Check agg.loom file size and the cluster-annotation table.
6. Decide LD-reference route (lift over vs GRCh38 reference).
7. Ask the teacher/mentor to confirm "no human subjects / no animal" answers on the form.
8. Identify who can serve as a scientific mentor for the genetics (the oncology mentor may not be
   the right person for statistical genetics).
