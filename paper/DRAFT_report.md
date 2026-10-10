# Which cells carry the genetic risk for ME/CFS, and is the pattern specific to the disease?

WORKING DRAFT. Methods and limitations describe what the code in this repository actually does, so they are
safe as a base. Introduction facts need their citations checked. Results are empty on purpose: nothing goes in
until the numbers exist. Rewrite everything in your own words before submitting; the prize needs your own writing
and you need to be able to explain every sentence.

Word budget to aim for: abstract under 300 words; report about 3,000 to 4,000 words plus figures.

---------------------------------------------------------------------------------------------------

## Abstract (write last, under 300 words)
Structure to follow: one sentence on why ME/CFS matters and is neglected; the gap (what is not known); what I
did (data, method, the new part); what I found (numbers from results/); what it means and its limits.
Do not use "first", "never", "groundbreaking", "cure" or "treatment".

---------------------------------------------------------------------------------------------------

## 1. Introduction
Points to cover, each with a checked citation:
1. ME/CFS: long-lasting illness with post-exertional malaise; WHO classifies it under neurological diseases
   (ICD-10 G93.3); mostly women (about 80% of DecodeME cases were female); often starts after an infection.
   Source for the sex ratio and onset: DecodeME preprint (Initial findings from the DecodeME GWAS, medRxiv 2025).
2. Burden: estimates of people affected (check the Institute of Medicine 2015 report and a recent review before
   quoting a number) and how little research money it gets compared with burden.
3. DecodeME: 15,579 cases and 259,909 population controls of European ancestry; eight genome-wide significant
   loci; SNP heritability on the liability scale about 0.095. The authors report enrichment in brain tissues.
4. The question this project asks: genetic signal for a trait is spread over many genes. Which cell types express
   those genes more than expected? Brain tissue is made of many cell types, so tissue-level enrichment cannot say
   which ones. A second problem: almost every polygenic brain trait lights up neurons, so finding neurons is weak
   evidence on its own. The useful question is whether ME/CFS differs from other traits.
5. What others have already shown, stated plainly: the DecodeME team and an informal forum analysis point to
   neurons. This project treats that as something to confirm, then tests specificity.

## 2. Methods

### 2.1 Data
- ME/CFS: DecodeME summary statistics (OSF project rgqs3), REGENIE output on GRCh38. Primary analysis `gwas_1`
  (15,579 cases, 259,909 UK Biobank controls). Robustness: `gwas_2` (same cases, a different control subset).
  Subtype and sex contrasts: infectious-onset, non-infectious-onset, female-only, male-only files.
- Only variants on the DecodeME quality-control list were used (MAF >= 0.01, INFO >= 0.9 and a DENTIST-style
  test), as their README recommends. An earlier run without this filter produced artefact signals in the
  pericentromeric regions of chromosomes 15 and 21 and was discarded (preregistration/DEVIATIONS.md).
- Cell types: Human Brain Cell Atlas v1.0 (Siletti et al. 2023), about 3 million nuclei from over 100 dissections of
  three donors, aggregated to 461 clusters in 31 superclusters.
- Comparison panel: [N] other GWAS (list and sources in data/panel.tsv), including positive controls
  (PGC3 schizophrenia, Alzheimer's disease) and a non-neural negative control (height).
- Reference genotypes: 1000 Genomes Phase 3 European individuals (n = 503).

### 2.2 Pre-registration
Gates, models, thresholds and the definition of an "ME/CFS-specific" cluster were written down and tagged in
git (`prereg-v1`, extension `prereg-v1.1`) before any ME/CFS cell-type result was opened. Changes after that are
logged with a date and reason in preregistration/DEVIATIONS.md.

### 2.3 Heritability gate
SNP heritability was estimated by LD score regression using the European LD scores (Bulik-Sullivan et al. 2015). The
standard software would not install on the Windows computer used (it depends on pysam), so the regression was
re-implemented in Python (scripts/ldsc.py), tested on simulated data, and checked against a published value: for
ME/CFS the observed-scale estimate converts to a liability-scale h2 close to DecodeME's reported 0.095 at a
prevalence of 0.5%. A GWAS was analysed further only if its heritability z-score was at least 4.

### 2.4 From SNPs to genes to cell types
DecodeME variants were mapped from GRCh38 to rsIDs by lifting positions to GRCh37 and matching position and allele
pair to the 1000 Genomes reference (97.6% matched). SNP p-values were combined into gene-level association
statistics with MAGMA 1.10 (SNPs within 35 kb upstream and 10 kb downstream of a gene; European LD reference).
For each cluster, a gene-level specificity score was computed: expression in counts per million, divided by the
gene's total across all 461 clusters, then rank-normalised within each cluster. 14,049 genes were kept (mean CPM at
least 1). A MAGMA gene-property regression asked whether genes more specific to a cluster are more strongly
associated with the trait.
- Model A conditions on the gene's mean expression across all clusters.
- Model B also conditions on its mean expression across neuronal clusters.
- Model C also conditions on the gene-level association of depression, BMI and insomnia.
Significance for model A used the Benjamini-Hochberg false discovery rate across 461 clusters (FDR < 0.05).

### 2.5 Specificity calibration
The same pipeline was run on every panel trait. For each trait, the 461 enrichment z-values were standardised
robustly (median and MAD across clusters). For ME/CFS in each cluster, the specificity score s is its standardised
value minus the panel mean for that cluster, divided by the panel standard deviation. A cluster is called
ME/CFS-specific only if FDR < 0.05 in model A, p < 0.05 in model B, and s >= 2.0.

### 2.6 Controls
Before interpreting ME/CFS, the pipeline had to reproduce known biology: schizophrenia enriched in excitatory
cortical neurons, Alzheimer's disease in microglia, and the three strongest height clusters being non-neuronal.

### 2.7 Translational ranking
Genes were ranked by the gene-level z multiplied by cluster specificity in the clusters enriched for the trait. The
ranking was benchmarked on migraine, rheumatoid arthritis and inflammatory bowel disease: do targets of approved,
selective drugs (Open Targets) rank higher than random genes matched on size and expression (10,000 permutations)?
Candidates for ME/CFS are labelled exploratory unless the benchmark passes. Everything is a hypothesis for
researchers, not a treatment suggestion.

## 3. Results
Facts and numbers only, taken from `results/` on 2026-10-10. Rewrite in your own words; keep every number traceable to the file named.

**3.1 Pipeline controls** (`results/G3_check.txt`, `results/method_benchmark.md`). Before any ME/CFS result was read, three control traits had to
reproduce known biology. Schizophrenia had an excitatory-neuron cluster at FDR < 0.05 (280 of 461 clusters significant in model A). Alzheimer's
disease was enriched in microglia (FDR 3e-7). Height's three strongest clusters were all fibroblasts. All three rules passed. Height also showed that
the continuous MAGMA test gives some false-positive neuronal clusters (13), so it is not used on its own to claim neuron specificity.

**3.2 Which cell types are enriched** (`results/clusters_decodeme_gwas_1.tsv`, `results/EVIDENCE_LEDGER.md`). 11 of 461 clusters pass FDR < 0.05 in
model A: deep-layer intratelencephalic (6), upper-layer intratelencephalic (2), amygdala excitatory (2) and one splatter cluster. The strongest
z is 3.58 (DLIT_136). All 11 meet the "confirmed" rule: S-LDSC, a gene-label permutation null (empirical p <= 0.002 at 1,000 permutations;
see the 5,000 run in `results/permutation_decodeme_gwas_1.tsv`), locus and leave-one-chromosome-out drops, and the second DecodeME analysis.
Medium spiny neurons, reported in earlier work, are not significant here (best z 2.8).

**3.3 Specificity against 19 other traits** (`results/clusters_decodeme_gwas_1.tsv`, column s). Six of the 11 clusters have s >= 2:
Amex_153 (s 2.98), Amex_175 (2.50), ULIT_121 (2.35), Splat_402 (2.41), DLIT_152 (2.03), DLIT_150 (2.01). The other five are enriched but shared with other traits.

**3.4 Conditioning** (`results/EVIDENCE_LEDGER.md`, Fig. 4). With mean neuronal expression also held constant (model B) no cluster passes FDR < 0.05
(best 0.067; all 11 have nominal p < 0.003). With depression, BMI and insomnia gene z-scores added (model C) 7 clusters pass, including DLIT_152,
DLIT_150, ULIT_121 and Amex_153.

**3.5 Robustness of the specificity** (`results/panel_robustness.tsv`, preregistration v1.8). Rule fixed beforehand: s >= 2 in at least 90% of
leave-one-trait-out panels and at least 80% of random half panels. Only the two amygdala excitatory clusters pass. The cortical clusters sit at
the cutoff (DLIT_150 and DLIT_152 have s of about 1.94 to 2.0 when one trait is dropped).

**3.6 Subtypes and sex** (`results/clusters_decodeme_gwas_1_*.tsv`; overlapping samples, descriptive). Infection-triggered onset: 55 clusters at FDR < 0.05;
non-infectious onset: 11; female-only 0; male-only 0. These counts follow sample size and power and are not tested as a difference.

**3.7 Where in the brain, and other layers** (Fig. 7, `results/layers_summary.txt`). The literature-region imaging convergence test is not significant
(10th percentile of the panel, permutation p 0.91). No Reactome pathway reaches FDR < 0.05. Of 154 whole-body cell types, 12 are at FDR < 0.05, including
non-neuronal ones, so the signal is not shown to be brain-only. Forward selection keeps one cluster for ME/CFS (two for schizophrenia).

**3.8 Candidate genes** (`results/candidate_genes_decodeme_gwas_1.tsv`, `results/benchmark.tsv`). The ranking placed known drug targets high for IBD
(median rank percentile 91%, permutation p 0.0002) but not for migraine or rheumatoid arthritis, so it fails the pre-registered benchmark (needs 3 of 3)
and the gene list (ISL1, CACNA1E, STT3B, DCC, PCDH17) is exploratory only.

## 4. Interpretation
Notes for the author to rewrite; the claims are limited to what the numbers support.
- **Question answered.** Earlier work says ME/CFS genetic risk sits in neurons. The question here was narrower: which cell types carry signal that
  other brain traits do not show? Across 461 clusters, 11 are enriched; most of that is a neuronal pattern shared with other traits.
- **What stands out.** Two amygdala excitatory clusters (Amex_153, Amex_175) keep a specificity score of 2 or more under every change to the
  comparison panel, and the permutation null (5,000 permutations) is significant for all 11 clusters. The cortical clusters sit at the cutoff.
- **What weakens it.** With neuronal expression also held constant (model B) no cluster passes FDR < 0.05, and conditioning on the principal
  components of 19 other traits leaves none at FDR < 0.05 (best 0.073). In bulk GTEx tissue the amygdala is not special among brain regions
  (bulk tissue can hide a small cell population). The second DecodeME analysis shares all cases, so it is not an independent replication.
- **What it means.** A hypothesis: amygdala excitatory neurons may be a cell population where ME/CFS genetics differs from other brain traits.
  It fits imaging reports of amygdala changes, but genetics here does not show mechanism, direction of effect or anything about treatment.
- **What is new.** The cell types were named informally before (forum analyses, Maccallini 2026). To my knowledge a panel-calibrated specificity
  test with robustness checks was not found in searches up to 10 October 2026 (to be repeated 17 October).
- **What would change my mind.** A larger independent ME/CFS GWAS where the amygdala clusters lose their specificity, or a single-cell
  amygdala atlas with more donors that places the signal elsewhere.

## 5. Pitfalls and limitations (true now, keep and extend)
- Cell-type enrichment is a statistical association, not proof that those cells cause the illness.
- The atlas comes from three post-mortem donors and covers the central nervous system; peripheral nerves, immune
  cells and autonomic ganglia, which matter in ME/CFS, are missing or thin.
- European ancestry only; 84% of cases are female, so male results are underpowered.
- Controls are UK Biobank participants (healthy-volunteer bias). Several panel GWAS also use UK Biobank, so
  controls overlap; this does not bias the per-trait enrichment but it matters for any cross-trait comparison.
- DecodeME case status comes from a questionnaire; cases were volunteers.
- GWAS-1 and GWAS-2 share all cases, so the second analysis is a robustness check, not an independent replication.
  The same holds for the female and infection-subset files.
- Specificity depends on the comparison panel; leave-one-out and alternative panels are reported as sensitivity.
- Gene prioritisation here uses association and cell-type specificity; it does not show direction of effect or causal
  mechanism. A genetic link to a gene does not mean a drug acting on it would help or be safe.
- The heritability code is a re-implementation (tested against simulations and one published figure).
- Some prior work could be missing from my literature searches, which stopped on [date].

## 6. Future work
Replicate in a larger ME/CFS cohort when summary statistics are released; add peripheral nervous system and
immune single-cell references; test top candidate genes in cell models with ME/CFS researchers.

## 7. Acknowledgements and data availability
Data: DecodeME (University of Edinburgh), Human Brain Cell Atlas, PGC, EBI GWAS Catalog, FinnGen, GTEx,
gnomAD, Open Targets. Code and pre-registration: github.com/palashraks-afk/AAN. List everyone who helped.

## References (to verify and format)
Bulik-Sullivan 2015 Nat Genet; de Leeuw 2015 PLoS Comput Biol (MAGMA); Siletti 2023 Science; DecodeME medRxiv 2025;
Minikel 2024 Nature; Watanabe 2019 Nat Commun (cell-type enrichment practice); Bryois 2020 Nat Genet (specificity
and MAGMA); Trubetskoy 2022 Nature (PGC3); Bellenguez 2022 Nat Genet; Yengo 2022 Nature; Kurki 2023 Nature (FinnGen).
