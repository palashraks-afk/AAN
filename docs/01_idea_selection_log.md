# Idea selection log (2026-10-07)

Brief from Palash: an idea that is (1) super advanced, (2) groundbreaking, (3) well-founded,
(4) new and novel, (5) genuinely unique. Later: must affect a human disease/condition with >50,000
people, have potential to inspire treatment, rate 4–5 on everything, human-centred, not necessarily
his field. Rules adopted: no invented ratings; "never done before" is only claimed after searching,
and even then as "not found".

## Ideas considered

| # | Idea | Verdict | Why |
|---|---|---|---|
| 1 | Data-leakage audit of EEG / MRI deep-learning classifiers | DROPPED | Saturated: Brookshire 2024 (EEG), Wen 2019 (ADNI MRI), 2025 schizophrenia-EEG preprint, July 2026 EEG-foundation-model stress test |
| 2 | NHANES cognition + blood markers / ML | DROPPED | Many NHANES cognition ML papers exist |
| 3 | NHANES smell test, cognition, mortality | DROPPED | Choi 2020 and others |
| 4 | EEG aperiodic exponent in Parkinson's/Alzheimer's: is it muscle? | SHELVED | Real gap (no shared-pipeline, muscle/fit-range test found), but moderate novelty; 2026 PD meta-analysis already pools open datasets; ds004504 is pre-cleaned to 0.5–45 Hz so cannot test EMG |
| 5 | Whole-fly-brain (FlyWire) in-silico disinhibition / seizure threshold | REJECTED BY USER ("no like") | No prior seizure study on the adult Shiu model found; closest = Huang 2018 (20k-neuron model, seizure-like firing only when wiring shuffled) and a 2026 inhibitory-hub preprint |
| 6 | Spaceflight mouse-brain transcriptome: cross-mission replication (NASA OSDR) | KEPT AS FALLBACK | RR-3 paper found only 11 shared genes vs 11 other datasets and called for replication; no pre-registered cross-mission brain replication found. Weaker human link. Datasets: OSD-33, 536, 352, 682, 685, 698, 699, 738, 751 |
| 7 | Parkinson's genetics in India (2025 GWAS, 4,806 cases / 6,364 controls) downstream analysis | BLOCKED | Zenodo record 14726796 files are RESTRICTED (access request form, approver/timeline unknown). No downstream cell-type/coloc/fine-mapping found |
| 8 | AI variant-effect predictors on gain- vs loss-of-function in epilepsy channels | DROPPED | Crowded: LoGoFunc 2023, GABA-A 2025, PIEZO1 2025, REVEL GoF-bias 2025, 2025 "distance metric" preprint |
| 9 | Neurology Neglect Atlas (GBD burden vs publications / NIH funding / trial locations) | KEPT AS STRETCH | No all-condition publication+trial-vs-burden study found; but NIH-funding-vs-burden work exists (PMID 41286576; Nature s44401-025-00048-x). Novelty 4, impact 3 |
| 10 | Snakebite neurotoxin antivenom design | NOT PURSUED | Search was refused by a safety filter; not worked around |
| 11 | Normal pressure hydrocephalus / essential tremor / IIH genetics | DROPPED | NPH FinnGen GWAS exists; ET GWAS (16,480 cases) + 2026 multi-omics exists; IIH cohorts tiny |
| 12 | **ME/CFS nervous-system cell-type map, specificity-calibrated + target dossier** | **SELECTED** | See PLAN.md |

## Why ME/CFS was selected

- Size: >1 million US; WHO-neurological; heavily neglected.
- Data: open summary statistics (no application) — verified file list and MD5 on OSF.
- Gap: brain-TISSUE enrichment published (13 brain tissues, p<0.05/54); an informal forum
  (S4ME) MAGMA run on Siletti clusters reported 24 significant neuronal cell types; I found NO
  specificity-calibrated, conditioned, subtype-aware, or translational cell-type analysis.
- Fit: the "measure, don't assume / negative control" style already used in Oncovision.

## Angles added to raise novelty/impact (each checked for prior art)

1. Specificity layer against a ≥20-trait panel (main novelty).
2. Conditioning on depression / BMI / insomnia shared genetics.
3. Autonomic/hypothalamic/brainstem resolution (whole-brain atlas, not only cortex).
4. Subtype contrasts: infection-triggered vs not; female vs all (DecodeME released these GWAS).
5. Pre-registered kill criteria; null result is publishable.
6. Translational dossier: direction-aware MR + druggability + safety + **benchmark on known drugs**.
7. Optional neglect-weighted panel of other neglected conditions.

## Honest limits of the search

~45 targeted searches over several hours. Absence of evidence is not proof. Prior-art overlaps found
for components (see docs/02). One search was blocked by a safety filter and dropped.
