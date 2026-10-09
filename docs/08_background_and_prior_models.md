# Background statistics, earlier research models, and how they connect to this project

Gathered 2026-10-08. Every number has a source; numbers I computed myself say so. Re-check the sources before quoting.

## A. The problem in numbers

### How many people, and how many are undiagnosed
| Measure | Value | Source |
|---|---|---|
| US cases, 2015 Institute of Medicine report | 836,000 to 2.5 million | IOM 2015 (via NCBI Bookshelf NBK284897, MEpedia summary) |
| US cases, updated for population growth | at least 1.5 million at the low end | Jason and Mirin, 2021 update (solvecfs.org PDF) |
| US cases, CDC | up to 3.3 million | cdc.gov/me-cfs/about |
| US adults, 2021-2022 national survey | 1.3% | NCHS figure reported via MEpedia (check original) |
| Prevalence in reviews | 0.42% (US study); 0.89% (13-country studies) | summarised in a long COVID review (PMC10654790 area); check original |
| Not diagnosed by a doctor | more than 9 in 10 (CDC); 84 to 91% (Medscape) | cdc.gov; Medscape 235980 |
| Share of cases that are female | about 80%; 84% in DecodeME | DecodeME preprint |
| ME/CFS cases among long COVID patients | 13 to 49% across studies; 58% in one cohort of 465; 8.4% in a Japanese clinic | Jason and Dorri (PMC9844405); PLOS ONE 0315385 |

### Cost
| Measure | Value | Source |
|---|---|---|
| Annual US economic cost, IOM 2015 | $17 to $24 billion | IOM 2015 |
| Updated for growth and inflation | $36 to $51 billion | Jason and Mirin 2021 |
| CDC current page | $18 to $51 billion in medical costs and lost income | cdc.gov |
One advocacy fact sheet gives far larger figures ($149 to $362 billion); I could not verify it, so it is not used.

### Research effort relative to burden (measured here and from sources)
| Measure | ME/CFS | Multiple sclerosis | Source |
|---|---|---|---|
| NIH research funding, about 2017 to 2019 | $15 million | $111 million | claims-data paper PMC6331450; ME/CFS Science 2019 comparison |
| US patients assumed | 1.73 to 3.75 million | 486,000 | same |
| NIH dollars per patient per year (my arithmetic) | $4 to $9 | $228 | `scripts/background_stats.py` |
| Gap per patient | about 26 to 57 times lower | | my arithmetic |
| PubMed papers 2021 to 2025 (title/abstract) | 1,833 | 30,173 | PubMed E-utilities, `results/background_pubmed_counts.tsv` |
| Papers per 1,000 patients | 0.49 to 1.06 | 62 | my arithmetic |
| Papers in 2025 | 402 | 6,073 | PubMed |
| Funding needed to match burden | about $203 million a year | | Mirin, Dimmock and Jason 2020 (Work 67) |
| Average NIH funding 1995 to 2014 | about $5 million a year, among the lowest relative to burden | | same |
Caveats: the funding and patient numbers come from different years and prevalence is uncertain, so the per-patient gap
is a range and an order of magnitude, not a precise figure. PubMed counts for "chronic fatigue syndrome" include some
non-ME/CFS fatigue papers. Funding has risen somewhat since; check current NIH RCDC values.

Figure: `figures/fig1_background.png` (made by `scripts/background_stats.py`).

### What is known about the biology (so the results can be set against it)
- PET with 11C-PK11195 (9 patients, 10 controls, 2014) found widespread brain neuroinflammation, including hippocampus,
  amygdala, thalamus, midbrain and pons (JNM 55:945). A later TSPO-PET/MRI abstract reported the putamen (ISMRM 2020).
  Newer TSPO work found no difference, and genotype-dependent ligand binding is a known problem.
- A 2023 meta-analysis of 65 neuroimaging studies found altered regions most often in cerebral cortex, with pooled
  hypoactivity in insula and thalamus (ScienceDirect S1568997223002185). A systematic review judged findings
  inconsistent (PLOS ONE 0232475).
- Autonomic evidence is the most replicated: POTS in 46% of patients against 19% of controls in one study (PMC9285104);
  meta-analysis of cardiac autonomic measures (PMC6824690). Brainstem and hypothalamic involvement is suggested by
  perfusion and MRI work and by orexin/HPA-axis reviews, mostly small studies.

## B. Earlier research models on ME/CFS genetics, compared

| Study or approach | Data and method | What it showed | What it did not do |
|---|---|---|---|
| DecodeME preprint (2025) | 15,579 cases, 259,909 controls; REGENIE GWAS; FUMA/MAGMA tissue enrichment; coloc with GTEx v10 | 8 loci; heritability 0.095; all 13 brain tissues enriched | no cell-type resolution, no test of specificity against other traits |
| DecodeME team talk (2026, not peer reviewed) | same data | heritability concentrates in neural tissue, neurons rather than glia | no published specificity or conditioning analysis seen |
| Forum MAGMA analysis (S4ME, informal) | DecodeME sumstats against 3,172 atlas region-cell combinations | 24 significant, all neurons; medium spiny neurons top (Phoenix Rising post) | no calibration, no conditioning, not peer reviewed |
| PrecisionLife combinatorial analysis (J Transl Med 2026) | combinations of SNPs in DecodeME cases | 259 core genes, over 40 candidate targets, links to long COVID | different statistical framework, no cell-type map |
| Gene-centric repurposing (PMC12940889) | gene sets from expression data and DecodeME MAGMA, compounds ranked by gene hits | candidate compounds | did not filter by direction of effect |
| Million Veteran Program with DecodeME meta-analysis | 19,470 cases, 699,111 controls | larger discovery sample | summary in the preprint's later analyses; no cell-type work found |
| Earlier GWAS and candidate-gene work | few thousand cases at most (for example 2,382 in UK Biobank 2022) | mostly underpowered | review: Hum Mol Genet 2020, "a critical review" |
| Neuroimaging (PET, MRI, fMRI) | small studies and meta-analyses | suggestive subcortical and brainstem findings | inconsistent across studies |
| This project | DecodeME sumstats, 461-cluster brain atlas, 23-trait panel, whole-body cell types, pathways | pre-registered calibration, conditioning, confirmation, benchmark, convergence | no individual-level data, no experiments, European ancestry only |

## C. Statistical methods for cell-type enrichment, compared

| Method | Idea | Strengths | Weaknesses | Used here |
|---|---|---|---|---|
| MAGMA gene-property (de Leeuw 2015) | regress gene-level association on a cell-type specificity score | fast, handles LD and gene size, can condition on other variables | gene-level step assumes a window model | primary test |
| LDSC-SEG, stratified LD score regression (Finucane 2018) | partition SNP heritability by annotations near specific genes | different statistical model; controls for baseline annotations | small annotations can inflate type 1 error (2021 medRxiv); full baseline model needed | simplified own version as confirmation (C1) |
| scDRS (Zhang 2022) | score each cell for disease relevance from GWAS genes | single-cell resolution, 74 traits | needs the single-cell object; different null | not run; cited |
| CELLECT / CELLEX, FUMA cell type (Watanabe 2019) | specificity scores plus MAGMA or LDSC | convenient web tools | choices of specificity metric matter | the specificity idea follows this line |
| Top-decile gene sets (Bryois 2020 style) | binary set of the most specific genes | simple, robust to outliers | loses information | planned sensitivity check (M1) |
Only the LDSC-SEG and scDRS rows were checked against search results this session; the other rows come from my own
knowledge of these papers and should be verified against the originals before citing.

### Method benchmark on this project's own control traits (M1)
The same positive and negative control traits are run through three methods: continuous MAGMA gene-property,
top-decile MAGMA set test, and the simplified S-LDSC. Reported: whether each recovers the expected cell classes
(schizophrenia, Alzheimer's, rheumatoid arthritis, height as non-neural), false-positive rate on height, agreement
between methods (Spearman across 461 clusters), and run time. Rules in `preregistration/PREREG_v1.4_method_benchmark.md`.

## D. How the background connects to the design

| Fact from the background | Design choice it explains |
|---|---|
| ME/CFS gets about 1/25 to 1/55 of the funding per patient that a comparable condition gets | open public data and secondary analysis; a low-cost way to add evidence |
| More than 9 in 10 patients undiagnosed; no objective biomarker | emphasis on biology that could guide biomarkers and targets, labelled as hypotheses |
| Case definitions differ between cohorts (Canadian Consensus in DecodeME; persistent tiredness in UK Biobank) | second GWAS, UK Biobank fatigue GWAS, subtype contrasts, and a stated limitation on definitions |
| 80 to 85% of cases are female | female-only analysis; male analysis passes the heritability gate only barely and is flagged underpowered |
| Illness often starts after an infection; large overlap with long COVID | infection-onset versus non-infection-onset contrast; comparison with long-COVID-related traits in the comorbidity map |
| Imaging findings are small and inconsistent | convergence test (E1) with a region list fixed before results, calibrated against other brain traits |
| Autonomic and brainstem evidence is the most replicated | brainstem, hypothalamus and sensory-ganglia cell types included in the atlas analysis; atlas limitation (no peripheral nerves) stated |
| Genetic support makes a drug target about 2.6 times more likely to succeed (Minikel 2024) | translational layer, benchmarked on diseases with known drugs, never presented as treatment advice |
| Prior ME/CFS work reports neurons but not specificity | the specificity calibration against a 23-trait panel is the main new test |

## E. More statistics on the problem (added 2026-10-08, evening)
| Measure | Value | Source and caveat |
|---|---|---|
| Time to diagnosis | median 12.0 years, mean 14.9 years among people diagnosed in 2022; median 1.5 years in 1987 | Solve ME/CFS registry analysis (solvecfs.org, 2023); weak trend, R2 = 0.39 |
| Full recovery in follow-up | median 5% (systematic review of 14 studies, Cairns and Hotopf 2005); 0 to 66% across 22 studies in a 2014 review | definitions and follow-up lengths differ; children recover more often (68% reported recovery by 10 years in one cohort) |
| Funding against burden | about 7% of the NIH funding its burden would predict, the lowest ratio in that analysis | Mirin, Dimmock and Jason 2020 (Work) |
| Projection with long COVID | 5 to 9 million US cases | projection, not a count (Jason and Mirin 2022 update) |

## F. How the other research approaches compare (immune, metabolomic, imaging, animal), added 2026-10-08
| Approach | What it can tell us | Typical sample | Main weakness | Where this project connects |
|---|---|---|---|---|
| GWAS (DecodeME) | inherited risk; genes come before illness | 15,579 cases, 259,909 controls | not the current disease state; UK Biobank controls; European ancestry | the data used here |
| Neuroimaging (PET, MRI, fMRI) | where the living brain differs | tens of patients (for example 26 cases, 18 controls) | small, inconsistent; brain-wide studies in other diseases need hundreds per group | independent check of the regions genetics points to (E1) |
| Immune profiling (cytokines, single-cell) | blood-cell differences | tens to hundreds (one multicentre study: 298) | no robust reproducible biomarker; sex and batch confound results | genetics tests whether signal sits in neural or immune cells (D2) |
| Metabolomics | chemical state of blood | 26 to 52 patients in several studies | no metabolite consistently altered; power weakened by many metabolites and few people | not used |
| Animal models | mechanism and treatment tests | lab cohorts | induced-inflammation models do not reproduce the human illness; protocols not standardised | later testing ground for candidate genes |
| This project | which cell types carry the inherited risk | same 15,579 cases | association only; three-donor atlas; CNS only | adds cellular resolution and checks it against the others |
Take-away: GWAS is the only ME/CFS approach with sample sizes in the tens of thousands; the others share small cohorts, mixed
case definitions and weak replication. A genetic signal that agrees with independent imaging and physiology is worth more
than any one alone, and one that disagrees is worth reporting. Sources: Solve ME/CFS 2023; Cairns and Hotopf 2005; reviews of
neuroimaging (PLoS One 2020), immune profiling, metabolomics and animal models (2025) found in searches of 2026-10-08;
verify before citing. No source compared replication rates head to head, so the ranking is qualitative.
