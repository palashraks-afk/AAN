"""Download the DecodeME summary statistics listed in data/manifest.tsv and verify them.

Not run automatically. Use --dry-run first. Downloads need Palash's approval.

    python scripts/download_decodeme.py --dry-run
    python scripts/download_decodeme.py --dest D:/AAN_data/decodeme
    python scripts/download_decodeme.py --only gwas_1.regenie.gz
"""
import argparse
import csv
import hashlib
import shutil
import sys
import time
from pathlib import Path

import requests

MANIFEST = Path(__file__).resolve().parent.parent / "data" / "manifest.tsv"
CHUNK = 1 << 20


def read_manifest():
    with open(MANIFEST, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url, out, expected_bytes, tries=8):
    # OSF answers 429 if you ask for several big files in a row, so back off and retry
    for attempt in range(tries):
        r = requests.get(url, stream=True, allow_redirects=True, timeout=60)
        if r.status_code == 429:
            wait = int(r.headers.get("Retry-After", 0)) or 30 * (attempt + 1)
            print(f"  429 from server, waiting {wait}s")
            r.close()
            time.sleep(wait)
            continue
        r.raise_for_status()
        with open(out, "wb") as fh:
            for block in r.iter_content(CHUNK):
                fh.write(block)
        r.close()
        break
    else:
        raise RuntimeError(f"gave up on {out.name} after {tries} tries")
    got = out.stat().st_size
    if got != expected_bytes:
        raise RuntimeError(f"{out.name}: got {got} bytes, expected {expected_bytes}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default="D:/AAN_data/decodeme")
    ap.add_argument("--only", nargs="*", help="file names to fetch (default: all)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rows = read_manifest()
    if args.only:
        rows = [r for r in rows if r["name"] in set(args.only)]
    dest = Path(args.dest)
    total = sum(int(r["bytes"]) for r in rows)
    probe = dest if dest.exists() else dest.parent if dest.parent.exists() else Path(dest.anchor)
    free = shutil.disk_usage(probe).free
    print(f"{len(rows)} files, {total / 1e9:.2f} GB; free on target drive: {free / 1e9:.1f} GB")
    for r in rows:
        print(f"  {r['name']:42s} {int(r['bytes']) / 1e6:9.1f} MB  {r['url']}")
    if args.dry_run:
        return 0
    if free < total * 1.2:
        print("Not enough free space; aborting.", file=sys.stderr)
        return 1

    dest.mkdir(parents=True, exist_ok=True)
    for r in rows:
        out = dest / r["name"]
        if out.exists() and out.stat().st_size == int(r["bytes"]) and md5_of(out) == r["md5"]:
            print(f"OK (already verified): {r['name']}")
            continue
        print(f"Downloading {r['name']} ...")
        fetch(r["url"], out, int(r["bytes"]))
        digest = md5_of(out)
        if digest != r["md5"]:
            raise RuntimeError(f"{r['name']}: md5 {digest} != {r['md5']}")
        print(f"Verified: {r['name']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
