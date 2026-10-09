"""Does the per-chromosome reference give the same gene-level results as the full reference? (chromosome 22 only)

Runs MAGMA's gene analysis for chromosome 22 of decodeme_gwas_1 twice, once against the full 22-million-SNP reference and
once against the split reference, and compares the gene z-scores.
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")
W = Path("D:/AAN_data/derived/magma/decodeme_gwas_1")
FULL = "C:/AAN_ref/g1000_eur"
SPLIT = "C:/AAN_ref/by_chr/g1000_eur_chr#CHR#"
annot = next(p for p in (W / "annot.genes.annot", W / "annot.genes.annot.txt") if p.exists())

for label, ref in (("full", FULL), ("split", SPLIT)):
    subprocess.run([str(MAGMA), "--bfile", ref, "--pval", str(W / "pval.txt"), "ncol=N", "--gene-annot", str(annot),
                    "--batch", "22", "chr", "--out", str(W / f"check_{label}")], check=True, stdout=subprocess.DEVNULL)

a = pd.read_csv(W / "check_full.batch22_chr.genes.out.txt", sep=r"\s+").set_index("GENE")
b = pd.read_csv(W / "check_split.batch22_chr.genes.out.txt", sep=r"\s+").set_index("GENE")
common = a.index.intersection(b.index)
print(f"genes: full {len(a)}, split {len(b)}, shared {len(common)}")
d = (a.loc[common, "ZSTAT"] - b.loc[common, "ZSTAT"]).abs()
print(f"max |dz| = {d.max():.2e}, correlation of z = {np.corrcoef(a.loc[common, 'ZSTAT'], b.loc[common, 'ZSTAT'])[0, 1]:.6f}")
print("same SNP counts:", bool((a.loc[common, "NSNPS"] == b.loc[common, "NSNPS"]).all()))
