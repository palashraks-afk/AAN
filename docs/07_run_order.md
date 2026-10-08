# Run order

Everything runs on Windows with Python 3.11 (virtualenv in `D:\AAN_data\venv`) and the MAGMA 1.10 Windows
binary. Data and large intermediates live in `D:\AAN_data`; the reference panel is split by chromosome
on the SSD in `C:\AAN_ref` because the D: drive is a slow hard disk.

| Step | Script | Output |
|---|---|---|
| 1 | `download_decodeme.py` | the six ME/CFS GWAS, variant list, INFO file (MD5 checked) |
| 2 | `build_rsid_map.py` | `decodeme_rsid_map.parquet` (GRCh38 variant id to rsID, 97.6% matched) |
| 3 | `build_specificity.py` | `gene_covar.txt`, `cluster_annotation.tsv`, `top10_sets.txt` |
| 4 | `split_reference.py` | per-chromosome copies of the 1000G EUR panel |
| 5 | `make_panel.py`, `download_panel.py` | comparison-trait GWAS (`data/panel.tsv`) |
| 6 | `ldsc.py` (tests in `tests/`) | SNP heritability and z-score for every GWAS, the G1 gate |
| 7 | `run_magma.py` | SNP to gene to models A and B for each trait |
| 8 | `check_controls.py` | G3: schizophrenia, Alzheimer's, height behave as expected |
| 9 | `calibrate.py` | `results/clusters_*.tsv`, specificity score and the T2 call |
| 10 | `make_conditioning_covar.py`, `gene_property.py --model C` | T3 conditioning |
| 11 | `get_known_targets.py`, `rank_targets.py` | benchmark X3 and the candidate list X1 |

## Things that went wrong and what fixed them (kept so the next run doesn't repeat them)
- The standard `ldsc` package needs `pysam`, which will not build on Windows. `scripts/ldsc.py` is a
  small re-implementation, checked by simulation and sanity-checked on schizophrenia (liability h2
  about 0.21 against a published 0.24).
- The original LDSC download link is dead. The European LD scores come from a Zenodo mirror (record
  8182036), MD5 matched.
- The GWAS Catalog legacy REST API is gone (HTTP 410). Study metadata comes from the downloadable study
  table and the FTP `harmonised_list.txt`.
- Open Targets removed `knownDrugs`; drugs now sit under `drugAndClinicalCandidates` with MONDO ids.
- A single MAGMA run on the full reference took hours. Splitting the reference by chromosome and using
  MAGMA's `--batch N chr` with `#CHR#` in the prefix cut a chromosome from about 8 minutes to 2.5.
- The Windows MAGMA writes the annotation as `annot.genes.annot.txt`.
- OSF rate-limits (HTTP 429) after a few large files; the downloader backs off and retries.
