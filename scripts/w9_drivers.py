"""PREREG v1.14 W9: does the amygdala signal rest on a few genes?"""
from pathlib import Path

import pandas as pd

import calibrate
import confirm_robustness as cr

DER = Path("D:/AAN_data/derived")
RES = Path("D:/AAN/results")
TARGET = "decodeme_gwas_1"
CLUSTERS = {"Amex_153": 153, "Amex_175": 175}


def main():
    wd = DER / "magma" / TARGET
    covar = pd.read_csv(DER / "gene_covar.txt", sep="\t")
    genes = pd.read_csv(wd / "genes.genes.out", sep="\t").set_index("GENE")["ZSTAT"]
    covar["zstat"] = covar["GENE"].map(genes)
    rows = []
    for name, cid in CLUSTERS.items():
        contrib = (covar["zstat"] * covar[f"c{cid}"]).dropna().sort_values(ascending=False)
        order = covar.loc[contrib.index, "GENE"].tolist()
        for k in (0, 1, 5, 20):
            keep = covar[~covar["GENE"].isin(order[:k])].drop(columns="zstat")
            path = wd / f"covar_drop_{name}_{k}.txt"
            keep.to_csv(path, sep="\t", index=False, float_format="%.5f")
            cr.run_model_a(wd, path, f"cluster_drop_{name}_{k}")

            df = pd.read_csv(wd / f"cluster_drop_{name}_{k}.gsa.out.txt", comment="#", sep=r"\s+").set_index("VARIABLE")
            p = float(df.loc[f"c{cid}", "P"])
            rows.append({"cluster": name, "genes_dropped": k, "p": p, "z": float(calibrate.z_from_p(p)),
                         "dropped": ";".join(str(g) for g in order[:k]) if k <= 5 else ""})
            path.unlink()
            print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "w9_drivers.tsv", sep="\t", index=False, float_format="%.4g")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
