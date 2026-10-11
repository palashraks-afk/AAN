"""PREREG v1.14 W2: false-positive rate of the T2 rule on null traits (gene rows permuted within size x expression bins)."""
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import calibrate
import confirm_robustness as cr

DER = Path("D:/AAN_data/derived")
RES = Path("D:/AAN/results")
TARGET = "decodeme_gwas_1"
N_NULL = 200


def main():
    rng = np.random.default_rng(20261010)
    wd = DER / "magma" / TARGET
    covar = pd.read_csv(DER / "gene_covar.txt", sep="\t")
    nsnps = pd.read_csv(wd / "genes.genes.out", sep="\t").set_index("GENE")["NSNPS"]
    covar["nsnps"] = covar["GENE"].map(nsnps)
    covar = covar.dropna(subset=["nsnps"]).reset_index(drop=True)
    covar["bin"] = (pd.qcut(covar["nsnps"], 5, labels=False, duplicates="drop") * 5
                    + pd.qcut(covar["avg_all"], 5, labels=False, duplicates="drop"))
    groups = [np.where(covar["bin"].to_numpy() == b)[0] for b in sorted(covar["bin"].unique())]
    cols = [c for c in covar.columns if c not in ("GENE", "nsnps", "bin")]

    names = calibrate.panel_names(calibrate.traits_with_results("A"))
    gate = calibrate.passes_heritability_gate(names)
    panel = [n for n in names if gate.get(n) is not False]
    prof = pd.DataFrame({t: calibrate.robust_standardise(calibrate.z_from_p(calibrate.load_model(t, "A")["P"]))
                         for t in panel})
    pm, psd = prof.mean(axis=1), prof.std(axis=1, ddof=1)

    rows = []
    for k in range(N_NULL):
        perm = np.arange(len(covar))
        for idx in groups:
            perm[idx] = rng.permutation(idx)
        # shuffle the 461 cluster columns together (genes keep their joint profile); avg_all stays with the gene
        shuf = covar[cols].iloc[perm].reset_index(drop=True)
        shuf["avg_all"] = covar["avg_all"].to_numpy()
        shuf.insert(0, "GENE", covar["GENE"].to_numpy())
        path = wd / "covar_null.txt"
        shuf.to_csv(path, sep="\t", index=False, float_format="%.4f")
        cr.run_model_a(wd, path, "null_run")
        df = pd.read_csv(wd / "null_run.gsa.out.txt", comment="#", sep=r"\s+")
        df = df[df["VARIABLE"].str.match(r"^c\d+$")]
        df["cid"] = df["VARIABLE"].str[1:].astype(int)
        df = df.set_index("cid").sort_index()
        z = calibrate.z_from_p(df["P"])
        zs = pd.Series(calibrate.robust_standardise(z), index=df.index)
        s = (zs - pm.reindex(df.index)) / psd.reindex(df.index)
        fdr = pd.Series(multipletests(df["P"], method="fdr_bh")[1], index=df.index)
        rows.append({"null": k, "n_fdr": int((fdr < 0.05).sum()), "n_specific": int(((fdr < 0.05) & (s >= 2)).sum()),
                     "n_p001": int((df["P"] < 0.001).sum()), "n_p001_s2": int(((df["P"] < 0.001) & (s >= 2)).sum())})
        if k % 10 == 9:
            print(k + 1, "nulls done", flush=True)
            pd.DataFrame(rows).to_csv(RES / "w2_null.tsv", sep="\t", index=False)
        path.unlink()
    out = pd.DataFrame(rows)
    out.to_csv(RES / "w2_null.tsv", sep="\t", index=False)
    print("share of nulls with >= 1 specific cluster (FDR<0.05 and s>=2):", (out.n_specific > 0).mean())
    print("share of nulls with any FDR<0.05 cluster:", (out.n_fdr > 0).mean())
    print("mean clusters per null with p<0.001 and s>=2:", out.n_p001_s2.mean(), "of", out.n_p001.mean(), "p<0.001")


if __name__ == "__main__":
    main()
