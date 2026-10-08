"""Build gene_covar_C.txt: the model-B covariate file plus gene-level z-scores of depression, BMI and insomnia.

The z-scores come from each trait's own MAGMA gene analysis (genes.genes.out, ZSTAT column), so the
three traits must have been run through run_magma.py first.
"""
from pathlib import Path

import pandas as pd

DERIVED = Path("D:/AAN_data/derived")
TRAITS = {"depression_broad": "z_depression_broad", "bmi": "z_bmi", "insomnia": "z_insomnia"}


def main():
    base = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t")
    for trait, col in TRAITS.items():
        genes = pd.read_csv(DERIVED / "magma" / trait / "genes.genes.out", sep=r"\s+")
        z = genes.set_index("GENE")["ZSTAT"]
        base[col] = base["GENE"].map(z)
        print(trait, "genes with z:", int(base[col].notna().sum()), "of", len(base))
    # genes MAGMA could not test for a trait get the median; max-miss in MAGMA drops a column above 5%
    base.to_csv(DERIVED / "gene_covar_C.txt", sep="\t", index=False, float_format="%.5f")
    print("wrote", DERIVED / "gene_covar_C.txt")


if __name__ == "__main__":
    main()
