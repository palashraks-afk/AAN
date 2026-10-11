# Questions a judge is likely to ask, with short honest answers (task M9, 2026-10-10)

Each answer is a starting point to rewrite in your own words and to be able to defend. The last column names the file that backs it. Numbers come from `results/` and `DISCOVERY.md`.

## The question and the data
| # | Question | Short honest answer | Backed by |
|---|---|---|---|
| 1 | Why ME/CFS, and why is it a neuroscience question? | WHO classes it as neurological (ICD-10 G93.3); over a million people in the US, mostly women, about 9 in 10 undiagnosed, and it gets a small fraction of research funding for its burden. The largest genetics study (DecodeME) found enrichment in brain tissue, so which brain cells is the next question. | docs/08, PLAN.md |
| 2 | Where does the data come from, and did you collect any? | Only public, de-identified summary statistics and public cell atlases. No patients, no samples, no animals, so the human-subjects answer is No. | docs/06 |
| 3 | What is a GWAS, and what does enrichment in a cell type mean? | A GWAS tests millions of DNA variants for association with a condition. Gene-level scores are then compared with how specifically each gene is expressed in a cell type; enrichment means risk genes are expressed more specifically in those cells than chance predicts. It is an association, not proof that those cells cause the illness. | docs/08 section C |
| 4 | How large is the study, and what are its weaknesses? | 15,579 cases and 259,909 controls, European ancestry, 84% female; cases are questionnaire-defined volunteers; controls are UK Biobank (healthy-volunteer bias). | docs/03, paper draft limitations |
| 5 | What do the heritability numbers say? | SNP heritability is modest (liability-scale 0.095 in the DecodeME paper; my own LD-score regression gave observed-scale 0.041, z 14), so genetics explains only part of the illness. | results/ldsc via ledger |

## The method
| # | Question | Short honest answer | Backed by |
|---|---|---|---|
| 6 | How did you link genes to cell types? | MAGMA gene-property analysis: gene-level association scores regressed on each cluster's expression specificity, conditioning on how broadly the gene is expressed; 461 clusters from the Human Brain Cell Atlas. | PREREG_v1, scripts/gene_property.py |
| 7 | Why is "neurons are enriched" not enough? | Almost every brain trait lights up neurons. So I asked whether ME/CFS exceeds what 19 other traits show for the same cluster (a specificity score s, threshold 2). | PREREG_v1 T2 |
| 8 | How do you know your pipeline works? | Three controls had to reproduce known biology before any ME/CFS result was read: schizophrenia in excitatory neurons, Alzheimer's in microglia (FDR 3e-7), height in non-neuronal cells (fibroblasts). It also passed a 4-of-5 check where each panel trait was the target in turn. | results/G3_check.txt, PREREG_v1.12 |
| 9 | Which of those checks did NOT work? | In the validation, schizophrenia's neuronal signal is shared with other brain traits, so none of its 280 clusters scores as specific; that is expected behaviour of a specificity score but I report it. | DISCOVERY.md (validation) |
| 10 | Why pre-register, and what changed? | Writing the plan first stops me choosing thresholds after seeing the data. Every change is dated in DEVIATIONS.md: the DecodeME files needed a quality-control filter, a p-value floor on the Mac, 1,000 versus 5,000 permutations, a panel-contamination bug I caught before reading the result. | DEVIATIONS.md, figures/fig_prereg_timeline.png |
| 11 | You re-implemented LD-score regression. Why trust it? | The standard package would not install; mine is tested on simulated data and gave a schizophrenia liability heritability near the published figure (0.21 vs 0.24). It is used for the heritability gate and a confirmation, not for the headline number. | tests/, docs/07 |
| 12 | Why not use a standard tool like FUMA? | Maccallini et al. used FUMA. I needed a panel of other traits and conditioning that FUMA's web interface does not do. | docs/11 |

## The result
| # | Question | Short honest answer | Backed by |
|---|---|---|---|
| 13 | What did you find? | 11 of 461 clusters are enriched (FDR < 0.05; excitatory cortical and amygdala neurons); six of them are more ME/CFS-specific than the other traits; two amygdala excitatory clusters stay specific however the comparison panel is changed. | EVIDENCE_LEDGER.md, panel_robustness.tsv |
| 14 | Is it a discovery? | No, not in the strong sense. The cell types were already pointed at informally, and the evidence is modest (FDR 0.03 to 0.05; z about 3 to 3.6). What is new in my searches is the panel-calibrated, pre-registered test. | docs/12 |
| 15 | Which of your own checks weakened it? | With neuronal expression held constant no cluster passes (best FDR 0.067); conditioning on 19 other traits' principal components, best FDR 0.073; bulk GTEx amygdala is not special; a second atlas and an independent fatigue cohort did not confirm it. | EXTENSION_LAYERS.md |
| 16 | Why did the independent checks not confirm it? | Each had limited power or a different tissue or phenotype: the second amygdala atlas has only 443 excitatory nuclei, the fatigue cohort is a broad symptom code, bulk tissue hides small cell groups. Nulls from such tests are weak evidence against, not refutations. | EXTENSION_LAYERS.md |
| 17 | What is the pain result? | The same six clusters show the same pattern in a large pain phenotype (mean z 2.46, p = 0.0001, all six positive), so "specific" holds against my 19 traits (which have no pain trait) but not against pain. It suggests shared circuitry for ME/CFS and chronic pain. | x_M4_related_conditions.tsv |
| 18 | Why the amygdala? Is that plausible? | The amygdala links stress, emotion and autonomic output, and imaging reports small, nominal amygdala differences in ME/CFS. My result does not show the amygdala causes anything; it says these cell types express risk genes unusually specifically. | docs/13, EXTENSION_LAYERS.md |
| 19 | What about the medium spiny neurons others reported? | In the main atlas they were not significant in my primary test (best z 2.8); in a second atlas indirect-pathway medium spiny neurons were nominal only. | DISCOVERY.md, x_ATLAS_tran.tsv |
| 20 | Did you find new drug targets? | No. A ranking of candidate genes was benchmarked on diseases with known drugs and passed only 1 of 3 (IBD), so it is labelled exploratory and I make no treatment claim. | results/benchmark.tsv |

## Limits and honesty
| # | Question | Short honest answer | Backed by |
|---|---|---|---|
| 21 | What would make your result wrong? | Association is not causation; the atlas is three post-mortem donors; peripheral nerves are barely covered; European ancestry; the second DecodeME analysis shares all cases. | paper draft section 5 |
| 22 | How many tests did you run, and is that p-hacking? | Many, which is why each one was pre-registered with a decision rule and why nulls are reported alongside the positive. The headline uses FDR across 461 clusters. The count and outcomes are in EXTENSION_LAYERS.md. | EXTENSION_LAYERS.md |
| 23 | Is the second DecodeME analysis a replication? | No. It uses the same cases with different controls; it tests stability, not independence. | DISCOVERY.md |
| 24 | Why are men and subtypes not analysed more? | Female-only and male-only analyses show no cluster at FDR < 0.05; the counts for infection-onset (55) versus non-infection (11) follow sample size, so no claim of a difference. | DISCOVERY.md F |
| 25 | What is your novelty claim, exactly? | "To my knowledge, a panel-calibrated specificity test of ME/CFS cell types was not found in searches of PubMed, Europe PMC (incl. medRxiv/bioRxiv), Google Scholar and citations on 10 Oct 2026." Repeat on 17 and 18 Oct. Not "first". | docs/11, docs/02 |
| 26 | Who else has done related work? | DecodeME (tissues), Maccallini et al. 2026 (cell types, FUMA), Lee et al. 2026 (fetal atlas), forum analyses (cell types), fibromyalgia Nature Medicine 2026 (neural cell types). | docs/13 |
| 27 | What did you do yourself, and what did tools do? | (Fill in your own words: you chose the question, designed the checks, ran and interpreted the analysis; list every person and tool that helped, accurately, in the application's "who helped" answer.) | TODO_PALASH.md |

## Looking forward
| # | Question | Short honest answer | Backed by |
|---|---|---|---|
| 28 | What would you do with more time or money? | Rerun on the larger meta-analyses when released (46,450 cases reported at a conference), add peripheral and autonomic single-cell data and spatial amygdala data, and test candidate genes in cell models with ME/CFS researchers. | paper draft section 6 |
| 29 | Could this help patients? | It can point researchers at cell types and shared circuitry (for example overlap with pain); it is not a treatment and should not change what anyone takes. | docs/06 |
| 30 | What would change your mind? | A larger cohort in which the amygdala clusters do not reach significance, or a calibrated analysis with a pain trait in the panel in which they stop being specific. | EXTENSION_LAYERS.md |
