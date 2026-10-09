# Setting up on another computer

Two ways: copy the finished data from the first machine (fast), or rebuild everything from the open sources (slow, about a day
of compute, but needs nothing from the first machine).

## Option A: copy the data
Copy `D:\AAN_data` (about 60 GB, external drive or cloud) and `C:\AAN_ref\by_chr` (about 6 GB) to the same relative layout, or
anywhere and replace the path prefixes (below). Then clone the repo, create the Python environment, and go to HANDOFF.md section 4.

## Option B: rebuild from scratch
1. Clone: `git clone https://github.com/palashraks-afk/AAN`
2. Python 3.11. Create a virtual environment and `pip install -r requirements.txt`; also `pip install openpyxl plotly nilearn nibabel pytest`.
3. Choose data folders and set them everywhere (the scripts hard-code `D:/AAN_data`, `C:/AAN_ref`, `D:/AAN`). Linux or macOS example:
   `grep -rl "D:/AAN_data" scripts | xargs sed -i 's#D:/AAN_data#/data/AAN_data#g'` and the same for `C:/AAN_ref` and `D:/AAN`.
   On Linux or macOS replace `magma.exe` by the MAGMA 1.10 binary for your system (below).
4. Tools (official sources):
   - MAGMA 1.10 and the gene location files (NCBI37.3 and NCBI38) and the 1000G European reference: https://cncr.nl/research/magma/
   - European LD scores (Zenodo mirror, MD5 e2f16343c4cfaa76caa7d0c03d26b489): https://zenodo.org/records/8182036
5. Data, in this order:
   - `python scripts/download_decodeme.py` (nine files, about 2.2 GB, MD5 checked; OSF rate-limits, the script retries)
   - atlas: aggregated loom link and `tables/cluster_annotation.xlsx` from https://github.com/linnarsson-lab/adult-human-brain
     (about 308 MB; do not download the 37 GB file)
   - `python scripts/make_panel.py` then `python scripts/download_panel.py` (about 8.5 GB) and PGC3 schizophrenia from figshare
     (article 19426775, European autosome file, MD5 6ebe2376f5cda972d37efa0f214c4df0)
   - FinnGen R13 files for G6_TRINEU, G6_CLUSTHEADACHE_WIDE, H8_MENIERE, G6_DYSTON from the public bucket
   - HPA single-cell-type table, Reactome pathways, GTEx v10 eQTL tar, gnomAD constraint (links in `docs/02` and `scripts/`)
6. Build: `python scripts/build_rsid_map.py`, `build_specificity.py`, `split_reference.py`, `build_discovery_covars.py`,
   `get_known_targets.py`, `python scripts/sldsc_lite.py --prepare`.
7. Run the heavy parts: `python scripts/run_overnight.py`, then `run_extension.py`, `run_sldsc_batch.py`.
8. Tests: `python -m pytest tests -q` (9 should pass).

## Keeping the machine awake
Long runs stop if the computer sleeps. Turn off sleep while they run (and keep it plugged in), or run them on a desktop.
