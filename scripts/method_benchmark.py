"""M1 (PREREG v1.4): do the three cell-type methods recover known biology on control traits, and do they agree?

Methods: (1) MAGMA continuous gene-property, model A; (2) MAGMA competitive set test on each cluster's top 10% most
specific genes; (3) the simplified S-LDSC. ME/CFS is not part of this benchmark.

    python scripts/method_benchmark.py
Writes results/method_benchmark.tsv and results/method_benchmark.md.
"""
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr
from statsmodels.stats.multitest import multipletests

import calibrate

DERIVED = Path("D:/AAN_data/derived")
RES = Path("D:/AAN/results")
MAGMA = Path("D:/AAN_data/tools/magma/magma.exe")
IMMUNE_TRAITS = {"rheumatoid_arthritis", "ibd"}


def prepare_inputs():
    sets = DERIVED / "top10_sets_magma.txt"
    if not sets.exists():
        with open(DERIVED / "top10_sets.txt") as src, open(sets, "w") as out:
            for line in src:
                name, genes = line.rstrip("\n").split("\t")
                out.write(name + " " + genes + "\n")
    covar = DERIVED / "gene_covar_avg.txt"
    if not covar.exists():
        pd.read_csv(DERIVED / "gene_covar.txt", sep="\t", usecols=["GENE", "avg_all"]).to_csv(covar, sep="\t", index=False)
    return sets, covar


def run_set_test(name, sets, covar):
    """Returns seconds taken (None if it was already done)."""
    workdir = DERIVED / "magma" / name
    if any((workdir / f"top10_sets.gsa.out{s}").exists() for s in (".txt", "")):
        return None
    t0 = time.time()
    subprocess.run([str(c) for c in [MAGMA, "--gene-results", workdir / "genes.genes.raw", "--set-annot", sets,
                                    "--gene-covar", covar, "--model", "condition-hide=avg_all", "analyse=sets",
                                    "--out", workdir / "top10_sets"]], check=True, stdout=subprocess.DEVNULL)
    return time.time() - t0


def load_gsa(name, prefix):
    for suffix in (".gsa.out.txt", ".gsa.out"):
        path = DERIVED / "magma" / name / f"{prefix}{suffix}"
        if path.exists():
            df = pd.read_csv(path, comment="#", sep=r"\s+")
            df = df[df["VARIABLE"].str.match(r"^c\d+$")]
            df["cluster_id"] = df["VARIABLE"].str[1:].astype(int)
            return df.set_index("cluster_id").sort_index()
    return None


def logged_seconds(name, prefix):
    log = DERIVED / "magma" / name / f"{prefix}.log"
    if not log.exists():
        return np.nan
    for line in log.read_text(errors="replace").splitlines():
        if "elapsed:" in line:
            h, m, s = line.split("elapsed:")[1].strip(" )").split(":")
            return int(h) * 3600 + int(m) * 60 + int(s)
    return np.nan


def sldsc_table(name):
    path = RES / f"sldsc_{name}.tsv"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    df["cluster_id"] = df["cluster"].str[1:].astype(int)
    return df.set_index("cluster_id").sort_index()


def expected_clusters(trait, ann):
    if trait == "schizophrenia":
        return ann.index[ann["supercluster"].isin(calibrate.EXCITATORY)]
    if trait == "alzheimer":
        return ann.index[ann["supercluster"] == "Microglia"]
    if trait in IMMUNE_TRAITS:
        return ann.index[ann["supercluster"] == "Miscellaneous"]
    return None


def main():
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    sets, covar = prepare_inputs()
    traits = [t for t in calibrate.traits_with_results("A")
              if not t.startswith("decodeme_") and t not in calibrate.NOT_IN_PANEL]
    rows, agree, seconds = [], [], {}
    for t in traits:
        seconds[t] = run_set_test(t, sets, covar)
        a, s, l = load_gsa(t, "cluster_A"), load_gsa(t, "top10_sets"), sldsc_table(t)
        z = {"magma_continuous": calibrate.z_from_p(a["P"]).reindex(ann.index)}
        p = {"magma_continuous": a["P"].reindex(ann.index)}
        if s is not None:
            z["magma_top10_set"] = calibrate.z_from_p(s["P"]).reindex(ann.index)
            p["magma_top10_set"] = s["P"].reindex(ann.index)
        if l is not None:
            z["sldsc"] = l["z"].reindex(ann.index)
            p["sldsc"] = l["p_one_sided"].reindex(ann.index)
        exp = expected_clusters(t, ann)
        for method, pv in p.items():
            fdr = pd.Series(multipletests(pv.fillna(1.0), method="fdr_bh")[1], index=pv.index)
            row = {"trait": t, "method": method, "n_fdr_sig": int((fdr < 0.05).sum()),
                   "n_neuronal_fdr_sig": int(((fdr < 0.05) & ann["neuronal"]).sum())}
            if exp is not None:
                row["recovers_expected"] = bool((fdr.loc[exp] < 0.05).any())
            if t == "height":
                top3 = z[method].sort_values(ascending=False).head(3).index
                row["top3_all_non_neuronal"] = bool((~ann.loc[top3, "neuronal"]).all())
            rows.append(row)
        if len(z) == 3:
            m = pd.DataFrame(z).dropna()
            agree.append({"trait": t, "cont_vs_set": spearmanr(m["magma_continuous"], m["magma_top10_set"])[0],
                          "cont_vs_sldsc": spearmanr(m["magma_continuous"], m["sldsc"])[0],
                          "set_vs_sldsc": spearmanr(m["magma_top10_set"], m["sldsc"])[0]})
    table = pd.DataFrame(rows)
    table["magma_continuous_seconds"] = table["trait"].map(lambda t: logged_seconds(t, "cluster_A"))
    table["top10_set_seconds"] = table["trait"].map(seconds)
    table.to_csv(RES / "method_benchmark.tsv", sep="\t", index=False)
    ag = pd.DataFrame(agree)
    ag.to_csv(RES / "method_benchmark_agreement.tsv", sep="\t", index=False, float_format="%.3f")

    controls = table[table["recovers_expected"].notna()]
    lines = ["# Method benchmark (M1, PREREG v1.4)", "",
             "Generated by scripts/method_benchmark.py. Control traits and panel only; ME/CFS is not part of it.", "",
             "## Do the methods recover the expected cell classes?", "",
             "| trait | method | recovers expected | FDR-significant clusters | neuronal FDR-significant |", "|---|---|:-:|---:|---:|"]
    for _, r in controls.iterrows():
        lines.append(f"| {r['trait']} | {r['method']} | {'yes' if r['recovers_expected'] else 'NO'} | "
                     f"{r['n_fdr_sig']} | {r['n_neuronal_fdr_sig']} |")
    h = table[table["trait"] == "height"]
    if len(h):
        lines += ["", "## Height (non-neural control): false positives", "",
                  "| method | neuronal clusters FDR < 0.05 | top three all non-neuronal |", "|---|---:|:-:|"]
        for _, r in h.iterrows():
            lines.append(f"| {r['method']} | {r['n_neuronal_fdr_sig']} | {'yes' if r.get('top3_all_non_neuronal') else 'no'} |")
    if len(ag):
        lines += ["", "## Agreement between methods (Spearman across 461 clusters)", "",
                  f"Median over {len(ag)} traits: continuous vs set test {ag['cont_vs_set'].median():.2f}, "
                  f"continuous vs S-LDSC {ag['cont_vs_sldsc'].median():.2f}, set test vs S-LDSC {ag['set_vs_sldsc'].median():.2f}."]
    lines += ["", "Rule (PREREG v1.4): a method that fails the schizophrenia or Alzheimer's control is not used to confirm ME/CFS."]
    (RES / "method_benchmark.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
