# Data

No data files are committed (`.gitignore` blocks them). Inputs live outside the repo in
`D:\AAN_data` because the C: drive has only ~13 GB free.

## Inputs

| Source | Files | Size | Access |
|---|---|---|---|
| DecodeME GWAS summary statistics (OSF `rgqs3`, folder "DecodeME Summary Statistics") | 9 files in `manifest.tsv` | ≈ 2.2 GB | Open listing; verify download works without login |
| Human Brain Cell Atlas v1.0 aggregated loom (github.com/linnarsson-lab/adult-human-brain) | `adult_human_20221007_Pool_Clean.agg.loom` + subcluster table | unknown — HEAD request first | Open (Dropbox / Google storage). DO NOT download the 37 GB full loom |
| Comparison-trait GWAS panel | ≥ 20 traits | TBD | Record URL, N, build, licence in `data/panel.tsv` |
| Positive-control drug-target lists | Open Targets / ChEMBL | small | Open API |

## Download policy
Downloads happen only after Palash approves them (state file, source, size). Run
`python scripts/download_decodeme.py --dry-run` first; it lists what it would fetch and checks free
space on the target drive. Every file is verified against `manifest.tsv` MD5 and bytes.

## Unverified
`README_shared_sum_stats.txt` could not be read (HTTP 429). Confirm REGENIE column names, allele
coding (ALLELE0/ALLELE1/A1FREQ), and sample-size column before writing any parser.
