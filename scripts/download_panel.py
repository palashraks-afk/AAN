"""Download the comparison-trait GWAS listed in data/panel.tsv (EBI GWAS Catalog harmonised files).

    python scripts/download_panel.py --dry-run
    python scripts/download_panel.py --skip fibromyalgia chronic_pain type1_diabetes

Each file is checked with a gzip integrity test after it lands. Already-complete files are skipped.
"""
import argparse
import gzip
import shutil
import sys
import time
from pathlib import Path

import pandas as pd
import requests

DEST = Path("D:/AAN_data/panel")
CHUNK = 1 << 20


def gzip_ok(path):
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(CHUNK):
                pass
        return True
    except (OSError, EOFError):
        return False


def fetch(url, out, tries=6):
    for attempt in range(tries):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                if r.status_code == 429:
                    time.sleep(30 * (attempt + 1))
                    continue
                r.raise_for_status()
                with open(out, "wb") as fh:
                    for block in r.iter_content(CHUNK):
                        fh.write(block)
            return True
        except requests.RequestException as err:
            print("  retry after error:", err)
            time.sleep(10 * (attempt + 1))
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip", nargs="*", default=[])
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    panel = pd.read_csv("D:/AAN/data/panel.tsv", sep="\t")
    panel = panel[~panel["label"].isin(args.skip)]
    if args.only:
        panel = panel[panel["label"].isin(args.only)]
    total = panel["bytes"].sum()
    print(f"{len(panel)} files, {total / 1e9:.2f} GB, free {shutil.disk_usage(DEST).free / 1e9:.0f} GB")
    if args.dry_run:
        print(panel[["label", "accession", "bytes"]].to_string(index=False))
        return 0

    DEST.mkdir(parents=True, exist_ok=True)
    for _, row in panel.iterrows():
        out = DEST / f"{row['label']}.h.tsv.gz"
        if out.exists() and out.stat().st_size == row["bytes"] and gzip_ok(out):
            print("have", row["label"])
            continue
        print("downloading", row["label"], f"{row['bytes'] / 1e6:.0f} MB", flush=True)
        if not fetch(row["url"], out) or not gzip_ok(out):
            print("FAILED", row["label"])
            continue
        print("ok", row["label"], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
