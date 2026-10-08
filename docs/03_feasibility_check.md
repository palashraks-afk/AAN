# Feasibility check — 2026-10-07

## Machine (Windows 10 Pro, 8 logical CPUs, ~16 GB RAM)
| Item | Result |
|---|---|
| Python | 3.11.9 |
| R | NOT installed (`Rscript` not found) |
| WSL | not reported (assume not available) |
| git | 2.53.0 |
| Python packages present | numpy 2.3.5, pandas 2.3.3, scipy 1.17.1, scikit-learn 1.9.0, matplotlib 3.10.8, seaborn 0.13.2, pyarrow 24, requests, pypdf, pdfplumber |
| Python packages MISSING | statsmodels, scanpy, anndata, h5py, pyliftover, polars, loompy |
| Genetics tools | MAGMA, LDSC, plink, coloc, TwoSampleMR: none installed |
| Disk | **C: 13 GB free (95% full)**; **D: 599 GB free** → everything goes on D: |
| GitHub CLI | `gh` not installed; plain `git` works; repo cloned to D:\AAN |

Implications: put data in `D:\AAN_data`; install MAGMA (standalone binary) and a Python-3 LDSC port;
install R + coloc/TwoSampleMR **or** implement ABF colocalisation in Python; a day of setup is planned.
RAM 16 GB is enough for summary statistics + the aggregated loom; NOT enough for the 37 GB full loom.

## DecodeME summary statistics (OSF node rgqs3, folder "DecodeME Summary Statistics")
Listed via the OSF API; downloads were NOT performed (permission not yet given for downloads).
Total ≈ 2.2 GB.

| File | Bytes | MD5 | Download |
|---|---:|---|---|
| gwas_1.regenie.gz | 318,442,179 | eabd3c06ffdeb2ec6382bfa67eed7f37 | https://osf.io/download/v4w8g/ |
| gwas_2.regenie.gz | 317,475,934 | 5c9338e5529a97433796c488bf6552e4 | https://osf.io/download/stp4h/ |
| gwas_1_female.regenie.gz | 450,200,982 | ca9b264e96d663dbd72afbe14fb0c35f | https://osf.io/download/3sqg4/ |
| gwas_1_male.regenie.gz | 305,246,888 | c3c3ea30311abb3619562c0510d56e7f | https://osf.io/download/z7dsq/ |
| gwas_1_infectious_onset.regenie.gz | 318,263,267 | 7429c0261d09cecc7034c98b7ed58c4c | https://osf.io/download/sz3fq/ |
| gwas_1_non_infectious_onset.regenie.gz | 316,049,034 | 68fb5ef38cd5771fec6903b881ac1230 | https://osf.io/download/zpjn8/ |
| gwas_qced.var.gz | 38,725,108 | b62e4dc634627b461c666bc444e9b0bf | https://osf.io/download/6uj5x/ |
| imputed.info.gz | 145,483,631 | 01943771c392ae2ba2c69d2b5f4a27bd | https://osf.io/download/h5edw/ |
| README_shared_sum_stats.txt | 1,914 | e1151ba8d0c40f53ca273de19071c9d5 | https://osf.io/download/axp4k/ |

Not yet verified: README contents (OSF returned HTTP 429 twice), exact REGENIE columns, whether
downloads need login (listing was anonymous; downloads redirect to files.de-1.osf.io).
Other OSF folders (not needed): Non-DecodeME-FMS-GWAS, UK Biobank - Samms and Ponting,
DecodeME Questionnaires; PDFs: GWAS analysis plan v1, questionnaire v6, regressionResults.pdf.
**GWAS-Analysis-Plan-v1.pdf is worth reading** (pre-registered plan from the DecodeME team).

## Facts extracted from the DecodeME preprint (read from the PDF)
- Cases: 26,901 completed the questionnaire (84% female); 21,620 met the study criteria (Canadian
  Consensus and/or required symptoms incl. post-exertional malaise); 18,051 samples received;
  15,579 cases after genotyping/QC and European-ancestry restriction.
- GWAS-1: 15,579 vs 259,909 (≈1:17). GWAS-Female 12,833 vs 218,949. GWAS-Male 2,746 vs 40,960.
  GWAS-2: all 15,579 cases vs 155,790 UKB controls not used in GWAS-1-style comparisons.
- 8,835,520 variants, genomic inflation λ = 1.003–1.066, GRCh38, REGENIE (Firth).
- h2 liability = 0.095 (SD 0.006). Eight genome-wide-significant loci overall.
- MAGMA: 13 genes of 18,637 significant; gene-tissue enrichment significant in all 13 brain tissues
  of 54 tissues tested.
- Controls are matched to the female:male ratio and ancestry of cases.

## Cell-type reference
- linnarsson-lab/adult-human-brain: full loom ~37 GB (do NOT download); aggregated
  `adult_human_20221007_Pool_Clean.agg.loom` (size not stated — HEAD request first);
  `Neurons.h5ad` / `Nonneurons.h5ad` (sizes not stated); subcluster annotation table in `tables`.
  Hosted on Dropbox and storage.googleapis.com/linnarsson-lab-human/. Code BSD-2.
- Cluster counts differ by source: 30/31/34 superclusters, 461 clusters, 3,313 subclusters. Fix
  ONE partition in the PREREG and report which.

## Other facts checked
- Zenodo 14726796 (Indian PD GWAS): files restricted (request form).
- FinnGen R13 release (Feb 2025): 2,755 endpoints, 500,186 participants; summary stats free from
  cloud storage; R14 planned Feb 2026 — check which is current.
- ADSP FunGen xQTL atlas (NIAGADS NG00184): open access, last release 2026-06-18, ~2 TB — too big.

## Risks found
1. Tool setup (MAGMA/LDSC/R) may eat 1–2 days → fallback to FUMA web server.
2. Build mismatch GRCh38 (DecodeME) vs GRCh37 (typical LD panels).
3. Comparison-trait GWAS availability unverified; many are UKB-based → control overlap.
4. Cluster-level power: 461 clusters × FDR; subcluster level (3,313) will be underpowered.
5. Time: 13 days to deadline; core path must be done by 13 Oct.
6. Novelty could be pre-empted by a paper not found.
