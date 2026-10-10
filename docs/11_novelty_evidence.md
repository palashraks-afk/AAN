# Novelty: the evidence behind the grade (3 now, 4 if the conditions are met)

Grade: **New and novel = 3 now; 4 only if the conditions in section 5 are met.** It was recorded as 4 on the morning of 2026-10-09
and lowered the same day, after the Mac session's re-search found (and this session then read in full) Maccallini et al. 2026
(section 2b). That preprint already reports cell-type enrichment for ME/CFS, a forward-selection test for independent signals and
a brain-versus-peripheral tissue result. A 5 is not available: the neuron result is a replication and no search proves that
nothing overlaps.
Searches were run 2026-10-07 to 2026-10-09 with a web search tool, not a library database. Repeat them in PubMed, bioRxiv and
medRxiv on **17 October** and record the result at the bottom of this file.

## 1. What the project claims is new (and what it does not)
| Layer | Claimed as new? | Pre-registration |
|---|---|---|
| ME/CFS signal is in neurons | **No, a replication** | PREREG v1 T1 |
| Is the neuronal signal specific to ME/CFS, tested against 23 other traits | **Yes, not found elsewhere** | v1 T2 |
| Same result after conditioning on generic neuronal expression and on depression, BMI, insomnia genetics | Yes | v1 models B and C |
| Two independent enrichment methods plus a gene-shuffling null plus locus and chromosome drop agree | Yes (a confirmation standard, not a biological finding) | v1.3 |
| How many independent cell-type axes there are (forward selection) | **Partly covered**: Maccallini 2026 ran FUMA forward selection and conditional analysis on ME/CFS. Ours adds the same test on schizophrenia, Alzheimer's and rheumatoid arthritis for comparison | v1.5 L1 |
| Direct region-level test and agreement with imaging literature, calibrated against other brain traits | Yes | v1.3 E1, v1.5 L2 |
| Brain versus whole-body cell types | **Mostly covered at tissue level** (Maccallini: GTEx, no peripheral tissue significant). Only the cell-type resolution (154 whole-body cell types) is new | v1.2 D2 |
| Pathways that still add after the cell type is accounted for | Partly: Maccallini tested 17,009 MSigDB sets (GO, Reactome) without conditioning on cell type; the conditioning step is new | v1.2 D3 |
| Infection-triggered versus other onset, at cell-type level | Yes (descriptive only, samples overlap) | v1 T5 |
| Cellular comorbidity map with long COVID and others | Yes | v1.2 D4, v1.6 |
| Cell-type-restricted target ranking validated on diseases with known drugs | Yes | v1 X1 to X3 |
| Neglected-condition extension | Exploratory; most conditions are expected to fail the heritability gate | v1.1 |

## 2. Closest prior work, layer by layer
| Prior work | What it did | What it did not do | Link |
|---|---|---|---|
| **Maccallini et al. 2026 (Research Square preprint, read in full)** | meta-analysis of DecodeME and the Million Veteran Program (19,470 cases, 699,111 controls); FUMA cell-type module (MAGMA gene-property, conditioning on average expression) on the Siletti, Seeker and DropViz brain atlases; FUMA forward selection for independent signals; 17,009 MSigDB gene sets; GTEx tissues (no peripheral tissue significant); eccentric medium spiny neurons (claustrum) and cerebellar white-matter glutamatergic neurons; glutamatergic synapses most specific | see section 2b | https://zenodo.org/records/20204356 and DOI 10.21203/rs.3.rs-9702020/v1 |
| DecodeME preprint (Aug 2025) | 15,579 cases GWAS, 8 loci, MAGMA tissue enrichment: all 13 brain tissues | no cell-type resolution, no specificity test | https://www.medrxiv.org/content/10.1101/2025.08.06.25333109v1.full |
| Ponting talk, International ME/CFS Conference 2026 | says heritability concentrates in neural tissue, neurons rather than glia; final paper still a preprint | not a peer-reviewed analysis; no calibration shown | https://events.mecfs-research.org/en/events/conference_2026/videos/chris-ponting-decodeme-genetics and https://www.meresearch.org.uk/highlights-from-the-international-me-cfs-conference-2026-including-decode-me-pet/ |
| Phoenix Rising forum post (21 Jul 2026) | informal cell-type run, medium spiny neuron as a top hit | informal; no calibration or conditioning | https://forums.phoenixrising.me/threads/cell-and-tissue-enrichment-in-me-cfsposted-on-july-21-2026.94698/ |
| S4ME forum MAGMA analysis | DecodeME against 3,172 atlas region-cell combinations, 24 significant, all neurons | no calibration, no conditioning, not peer reviewed | https://s4me.info/threads/using-magma-on-me-cfs-genetic-data.50465/ |
| PrecisionLife combinatorial analysis (J Transl Med 2026) | 259 core genes, over 40 candidate targets | different framework, no cell-type map | https://link.springer.com/article/10.1186/s12967-026-08167-1 |
| Gene-centric repurposing (PMC12940889) | compound ranking from DecodeME gene sets | no direction of effect, no cell-type restriction | https://pmc.ncbi.nlm.nih.gov/articles/PMC12940889/ |
| Long COVID GWAS (23andMe preprint) | genetic links to chronic fatigue, fibromyalgia, depression | no cell-type layer | https://www.medrxiv.org/content/10.1101/2024.10.07.24315052v1 |
| Larger ME/CFS meta-analysis (newsletter, unverified) | reported 46,450 cases, 2.46 million controls, 20 cohorts | I could not verify venue or content | https://connectivedigest.substack.com/p/summer-2026-mecfs-and-long-covid |
| scDRS (Zhang 2022) | single-cell disease relevance across 74 traits | no neglected neurological conditions, no ME/CFS | https://www.nature.com/articles/s41588-022-01167-z |
| Single-nucleus TWAS of 12 brain disorders (Nature 2026) | cross-disorder cell-type convergence | not ME/CFS or neglected conditions | https://www.nature.com/articles/s41586-026-10836-6 |
| Cell-type genetics of chronic pain (brain and DRG) | partitioned heritability by cell type | a different disease | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12700546/ |
| Tinnitus brain cell mapping (2025); restless legs (2024); essential tremor (2024) | cell-type work for those conditions | not ME/CFS; the reason those conditions are not in the extension | https://pmc.ncbi.nlm.nih.gov/articles/PMC12687797/ |
| ME/CFS imaging: PET neuroinflammation (2014), TSPO-PET putamen (2020 abstract), brainstem and hypothalamus work | independent regional evidence | small and inconsistent; used here as the convergence test, not as proof | https://jnm.snmjournals.org/content/55/6/945 and https://archive.ismrm.org/2020/3057.html |
| Minikel et al. 2024 | genetic support makes a drug target about 2.6 times more likely to succeed | not ME/CFS; justifies the translational layer | https://pmc.ncbi.nlm.nih.gov/articles/11096124/ |

## 2b. Maccallini et al. 2026: what it does and does not contain (full text read 2026-10-09)
Does: GWAS meta-analysis; MAGMA gene-set, GTEx tissue and FUMA cell-type analyses; step-2 forward selection and step-3 conditional
analysis across atlases (it concludes the many eMSN entries are one distributed signal); replication against a rare-variant gene
module (HEAL2). Cell types: eccentric medium spiny neurons (claustrum, S-47, Bonferroni P = 0.017) and cerebellar white-matter
glutamatergic neurons; 107 atlas datasets, 2,082 cell-type tests.
Does not (no match in the text; searched for: panel, other trait, depression, BMI, insomnia, negative control, permutation, LDSC,
partitioned, drug, repurpos, sex, infection): comparison of ME/CFS with other traits or any specificity test; conditioning on
depression, BMI or insomnia genetics; a second enrichment method; a permutation null; locus or chromosome drop; agreement with
imaging; sex or infection-onset analyses (their sources have no sex-specific statistics); benchmark-validated target ranking.
Their meta-analysis summary statistics are public (Zenodo, 312 MB). They include DecodeME cases, so they are not an independent
replication.

## 3. Searches run (terms and what came back)
| Date | Query (shortened) | Result |
|---|---|---|
| 10-07 | ME/CFS genetics cell type enrichment brain single-cell DecodeME MAGMA specificity neurons glia | no peer-reviewed analysis; forum posts and a conference talk |
| 10-07 | ME/CFS GWAS cell-type enrichment single-nucleus brainstem hypothalamus dorsal root ganglion autonomic Siletti | none |
| 10-07 | long COVID GWAS shared architecture ME/CFS cell types specificity vs depression, chronic pain, fibromyalgia | no cell-type result |
| 10-07 | DecodeME GWAS published journal 2026, final paper, cell types | still a preprint; neurons mentioned in a talk |
| 10-07 | ME/CFS drug repurposing genetic evidence DecodeME druggable genes Open Targets MR | gene-level repurposing only |
| 10-07 | Human Brain Cell Atlas Siletti GWAS cell-type enrichment | forum MAGMA run only |
| 10-08 | cross-disorder cell-type enrichment atlas neglected neurological conditions | snTWAS 2026, scDRS; none for neglected conditions |
| 10-08 | essential tremor / trigeminal neuralgia / cluster headache cell-type enrichment | none for the last two; essential tremor partial |
| 10-08 | tinnitus / Meniere / restless legs / narcolepsy cell-type enrichment | tinnitus and restless legs done; Meniere none |
| 10-08 | ME/CFS neuroimaging basal ganglia striatum hypothalamus brainstem | independent regional evidence found (used for E1) |
| 10-08 | ME/CFS autonomic brainstem hypothalamic evidence | autonomic evidence strongest |
| 10-08 | ME/CFS research approaches compared | no head-to-head comparison found |
About 45 targeted searches in total across the project, including earlier directions that were dropped.

## 4. Why 3 now, and why 4 is reachable
- Localisation, forward selection, the tissue-level brain-versus-body result and gene-set enrichment are all published (section 2b),
  so none of them can be claimed as new.
- Still not found anywhere: panel-calibrated specificity (is the signal specific to ME/CFS?), conditioning on depression, BMI and
  insomnia genetics, two-method confirmation plus permutation null and locus drop, a genetics-versus-imaging convergence test, sex
  and infection-onset contrast at cell-type level, benchmark-validated target ranking.
- The field is young (first large ME/CFS GWAS in 2025) and heavily under-researched.

## 5. Conditions for a 4, and what keeps it at 3
All of these must hold for a 4:
1. The control gate passes and the calibration layers actually run on the ME/CFS data.
2. They give a clear result: a specific, confirmed axis beyond what other brain traits show, or a clearly reported null.
3. The 17 October re-search finds nothing else that does the calibrated analysis.
It stays at 3 if any of these fails, or if a published version of Maccallini or the DecodeME paper adds cross-trait calibration.
Not available at all: a 5.

## 6. How to word it (use these, not stronger words)
- "ME/CFS is a young and under-researched field; the first large genome-wide study appeared in 2025."
- "Earlier work, including the DecodeME team's, points to neurons. I tested whether that neuronal signal is specific to ME/CFS or is
  shared by other polygenic brain traits, using a pre-registered comparison panel."
- "Cell-type enrichment in neurons, including striatal and cerebellar neurons, has been reported for ME/CFS (Maccallini et al. 2026;
  the DecodeME team). I tested whether that signal is specific to ME/CFS or shared by other polygenic brain traits."
- "To my knowledge this specificity test has not been reported for ME/CFS (searched PubMed, bioRxiv, medRxiv on [date])."
- Do not write: "first", "never done", "groundbreaking", "novel discovery" before the result exists.

## 7. 17 October re-search protocol
Databases: PubMed, bioRxiv, medRxiv, Google Scholar. Terms (each with "ME/CFS" or "myalgic encephalomyelitis" or "chronic fatigue
syndrome"): cell type, cell-type enrichment, single-nucleus, single-cell, MAGMA, LDSC, partitioned heritability, scDRS, specificity,
brainstem, hypothalamus, striatum, medium spiny, DecodeME, GWAS meta-analysis. Also: has the DecodeME paper been published, and
what cell-type analysis does it include? Record below.

### Log of re-searches
| Date | Searcher | Terms | New overlapping work found | Grade after |
|---|---|---|---|---|
| 2026-10-09 | Mac session found Maccallini; full text read here | see sections 2b and 3 | Maccallini et al. 2026 (cell types, independent signals, gene sets) | 3 (4 if section 5 holds) |
| 2026-10-17 | | | | |
