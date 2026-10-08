"""Split the 1000G EUR PLINK reference into one file set per chromosome.

MAGMA's chromosome batch mode reloads the whole 22-million-SNP reference for every chromosome, which
costs minutes per job. With one file set per chromosome (g1000_eur_chr1 ... g1000_eur_chr22) each job
only reads its own SNPs.

The .bed format is SNP-major: a 3-byte header, then ceil(n_individuals / 4) bytes per SNP, in the
same order as the .bim rows. So a chromosome is just a contiguous slice of both files.
"""
import shutil
from pathlib import Path

import pandas as pd

SRC = Path("C:/AAN_ref/g1000_eur")
OUT = Path("C:/AAN_ref/by_chr")
CHUNK = 64 * 1024 * 1024


def main():
    OUT.mkdir(exist_ok=True)
    n_ind = sum(1 for _ in open(f"{SRC}.fam"))
    row_bytes = (n_ind + 3) // 4
    bim = pd.read_csv(f"{SRC}.bim", sep="\t", header=None, dtype=str)
    with open(f"{SRC}.bed", "rb") as bed:
        header = bed.read(3)
        assert header == b"\x6c\x1b\x01", "not a SNP-major PLINK bed file"
        for chrom in [str(c) for c in range(1, 23)]:
            idx = bim.index[bim[0] == chrom]
            first, last = idx[0], idx[-1]
            assert len(idx) == last - first + 1, "chromosome rows are not contiguous"
            prefix = OUT / f"g1000_eur_chr{chrom}"
            bim.loc[first:last].to_csv(f"{prefix}.bim", sep="\t", header=False, index=False)
            shutil.copy(f"{SRC}.fam", f"{prefix}.fam")
            bed.seek(3 + first * row_bytes)
            remaining = len(idx) * row_bytes
            with open(f"{prefix}.bed", "wb") as out:
                out.write(header)
                while remaining:
                    block = bed.read(min(CHUNK, remaining))
                    out.write(block)
                    remaining -= len(block)
            print(f"chr{chrom}: {len(idx):,} SNPs", flush=True)


if __name__ == "__main__":
    main()
