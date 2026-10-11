# Reproducing this project (for someone who has never seen it)

## What it does in one paragraph
It asks which human brain cell types carry the inherited risk for ME/CFS, using the DecodeME GWAS summary statistics and a single-nucleus brain atlas, and whether the pattern is specific
to ME/CFS compared with 19 other traits. Every analysis was written down (git tags `prereg-v1` to `prereg-v1.16`, files in `preregistration/`) before it was run; changes are in `preregistration/DEVIATIONS.md`.

## Level 1: check the headline numbers in under a minute (no big data)
```bash
git clone https://github.com/palashraks-afk/AAN && cd AAN
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r requirements.lock        # exact versions used (requirements.txt has the unpinned list)
python -m pytest tests -q               # 14 checks: unit tests of the LD-score code and the committed result numbers
python scripts/x_figures.py             # rebuilds figures/fig_summary_panels.png and figures/fig_prereg_timeline.png from results/
```
Expected: all tests pass; the committed results show 11 clusters at FDR < 0.05 and six that are also specific (s >= 2).

## Level 2: rerun the full analysis (about a day on a laptop with an SSD, about 60 GB)
1. Follow `docs/09_new_machine_setup.md` (data sources, MD5 checks, MAGMA 1.10 binary for your system, European LD scores, 1000 Genomes reference).
2. Windows paths are written as `D:/AAN_data`, `C:/AAN_ref`, `D:/AAN`. On macOS or Linux run `scripts/sync_mac.sh`, which builds `~/AAN_mac` with your paths and the right MAGMA name, then run from there.
3. `python scripts/run_everything.py --dry-run` lists every step; `python scripts/run_everything.py` runs them in order, with the control gate (schizophrenia, Alzheimer's, height must reproduce known biology before any ME/CFS result is read).
4. Extension layers: `python scripts/x_pc_conditioning.py --validate --mecfs`, `x_gtex.py --build --run`, `x_atlas.py --build/--run`, `x_related.py`, `x_finngen_fatigue.py`, `x_m3_m5.py` (details in `EXTENSION_LAYERS.md`).

## Data availability statement (ready for a repository DOI)
All inputs are public: DecodeME summary statistics (OSF project rgqs3), Human Brain Cell Atlas v1.0 (Siletti et al. 2023), the comparison GWAS (GWAS Catalog accession numbers in `data/panel.tsv`; PGC3 schizophrenia on figshare),
FinnGen R13 public summary statistics, GTEx v8 median expression, the Tran et al. 2021 and Cao et al. 2020 atlases (CELLxGENE), HGNC gene groups and the MVP ME/CFS GWAS (GCST90479178).
No individual-level data or human participants were used. Code is MIT-licensed in this repository; result tables are in `results/`. Large inputs are not committed (see `data/README.md`).

## Known limits of the package
- Level 2 needs the large data files; they are not in the repository. - The own LD-score and S-LDSC code are re-implementations (tested against simulations and published numbers, see `tests/`).
- The numbers above were produced on two machines (Windows and macOS). The p-value floor of 1e-300 used for MAGMA input is in `DEVIATIONS.md`.
