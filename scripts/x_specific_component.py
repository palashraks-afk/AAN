"""PREREG v1.7 R1/R2: the ME/CFS-specific genetic component by gene-level residualisation.

    python scripts/x_specific_component.py --validate      # gate V (controls only, no ME/CFS)
    python scripts/x_specific_component.py --mecfs        # R1 and R2 for decodeme_gwas_1 and decodeme_gwas_2
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr
from statsmodels.stats.multitest import multipletests

import calibrate

DERIVED = Path("D:/AAN_data/derived")
RESULTS = Path("D:/AAN/results")
EXCLUDE = {"fibromyalgia", "chronic_pain", "type1_diabetes", "finngen_malaise_fatigue"}  # related conditions never sit in the null panel (PREREG v1.6, v1.11)
SEED = 20261010


def load_covar():
    cov = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t").set_index("GENE")
    spec_cols = [c for c in cov.columns if c.startswith("c") and c[1:].isdigit()]
    spec_cols = sorted(spec_cols, key=lambda c: int(c[1:]))
    return cov, spec_cols


def load_genes(name):
    f = DERIVED / "magma" / name / "genes.genes.out"
    g = pd.read_csv(f, sep=r"\s+").drop_duplicates("GENE").set_index("GENE")
    return g[["NSNPS", "NPARAM", "ZSTAT"]]


def panel_labels():
    """The pre-registered panel is exactly the labels of data/panel.tsv (schizophrenia is the PGC3 file); nothing analysed later can join it."""
    return set(pd.read_csv(Path(__file__).resolve().parent.parent / "data" / "panel.tsv", sep="\t")["label"])


def panel_traits():
    names = [n for n in calibrate.panel_names(calibrate.traits_with_results("A")) if n in panel_labels()]
    gate = calibrate.passes_heritability_gate(names)
    keep = [n for n in names if n not in EXCLUDE and gate.get(n) is not False]
    return keep


def kaiser_pcs(Z, k_cap=6, k=None):
    Zs = (Z - Z.mean(0)) / Z.std(0)
    u, s, vt = np.linalg.svd(Zs, full_matrices=False)
    ev = s ** 2 / (Zs.shape[0] - 1)
    kk = k if k is not None else int(min(k_cap, max(1, (ev >= 1).sum())))
    pcs = u[:, :kk] * s[:kk]
    return pcs, kk, ev[:kk].sum() / ev.sum()


def residualise(y, X):
    X1 = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return y - X1 @ beta, 1 - np.var(y - X1 @ beta) / np.var(y)


def bins_for(cov_df, n_size=10, n_expr=10):
    a = pd.qcut(cov_df["NSNPS"].rank(method="first"), n_size, labels=False)
    b = pd.qcut(cov_df["avg_all"].rank(method="first"), n_expr, labels=False)
    return (a * n_expr + b).to_numpy()


def cluster_stats(r, S_res, ss, bins, n_perm, rng):
    """Covariate-adjusted slope of r on every cluster's specificity, and a within-bin permutation null."""
    obs = S_res.T @ r / ss
    groups = [np.flatnonzero(bins == b) for b in np.unique(bins)]
    null = np.empty((n_perm, S_res.shape[1]))
    for i in range(n_perm):
        rp = r.copy()
        for idx in groups:
            rp[idx] = r[rng.permutation(idx)]
        null[i] = S_res.T @ rp / ss
    mu, sd = null.mean(0), null.std(0)
    z = (obs - mu) / sd
    p_emp = (1 + (null >= obs).sum(0)) / (n_perm + 1)
    return obs, z, p_emp


def prepare_design(cov, spec_cols, genes):
    X = cov.loc[genes, ["avg_all", "avg_neuron"]].to_numpy()
    X1 = np.column_stack([np.ones(len(genes)), X])
    S = cov.loc[genes, spec_cols].to_numpy()
    S_res = S - X1 @ np.linalg.lstsq(X1, S, rcond=None)[0]
    return S_res, (S_res ** 2).sum(0)


def run_target(target_z, target_meta, pcs, cov, spec_cols, genes, n_perm, rng):
    X = np.column_stack([pcs, np.log(target_meta["NSNPS"]), np.log(target_meta["NPARAM"])]) if pcs is not None \
        else np.column_stack([np.log(target_meta["NSNPS"]), np.log(target_meta["NPARAM"])])
    r, r2 = residualise(target_z, X)
    S_res, ss = prepare_design(cov, spec_cols, genes)
    bins = bins_for(pd.DataFrame({"NSNPS": target_meta["NSNPS"], "avg_all": cov.loc[genes, "avg_all"].to_numpy()}), 10, 10)
    obs, z, p = cluster_stats(r, S_res, ss, bins, n_perm, rng)
    fdr = multipletests(p, method="fdr_bh")[1]
    return pd.DataFrame({"cluster_id": np.arange(len(z)), "slope": obs, "z": z, "p_emp": p, "fdr": fdr}), r2


def common_genes(cov, traits):
    sets = [set(load_genes(t).index) for t in traits]
    g = set(cov.index)
    for s in sets:
        g &= s
    return sorted(g)


def raw_profile(name, cov, spec_cols, n_perm, rng):
    g = load_genes(name)
    genes = sorted(set(g.index) & set(cov.index))
    gg = g.loc[genes]
    return run_target(gg["ZSTAT"].to_numpy(), gg, None, cov, spec_cols, genes, n_perm, rng)[0]


def validate(n_perm):
    rng = np.random.default_rng(SEED)
    cov, spec_cols = load_covar()
    panel = panel_traits()
    out = []
    print("panel traits:", panel)
    ok_v1 = True
    for t in ("schizophrenia", "alzheimer", "height"):
        if t not in panel and t not in calibrate.traits_with_results("A"):
            print("missing", t)
            continue
        prof = raw_profile(t, cov, spec_cols, n_perm, rng)
        magma_z = calibrate.z_from_p(calibrate.load_model(t, "A")["P"].to_numpy())
        rho = spearmanr(prof["z"], magma_z)[0]
        print(f"V1 {t}: Spearman(OLS-permutation z, MAGMA z) = {rho:.3f}")
        out.append({"check": f"V1_{t}", "value": rho, "pass": rho >= 0.8})
        ok_v1 &= rho >= 0.8
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t")
    for t, kind in (("alzheimer", "Microglia"), ("height", "neuron")):
        others = [x for x in panel if x != t]
        genes = common_genes(cov, others + [t])
        Z = np.column_stack([load_genes(x).loc[genes, "ZSTAT"].to_numpy() for x in others])
        pcs, k, expl = kaiser_pcs(Z)
        meta = load_genes(t).loc[genes]
        res, r2 = run_target(meta["ZSTAT"].to_numpy(), meta, pcs, cov, spec_cols, genes, n_perm, rng)
        res = res.merge(ann, left_on="cluster_id", right_index=True, how="left") if "cluster_id" not in ann.columns \
            else res.merge(ann, on="cluster_id", how="left")
        sc = [c for c in res.columns if "supercluster" in c.lower()][0]
        if kind == "Microglia":
            hit = ((res[sc] == "Microglia") & (res["fdr"] < 0.05)).any()
            print(f"V2 alzheimer residual: microglia cluster at FDR<0.05 = {hit} (k={k}, PCs explain {r2:.2f} of trait variance)")
            out.append({"check": "V2_alzheimer_microglia", "value": float(hit), "pass": bool(hit)})
        else:
            neuro = res[~res[sc].isin(calibrate.NON_NEURONAL)]
            hit = (neuro["fdr"] < 0.05).any()
            print(f"V3 height residual: neuronal cluster at FDR<0.05 = {hit}")
            out.append({"check": "V3_height_no_neurons", "value": float(hit), "pass": not hit})
    df = pd.DataFrame(out)
    df.to_csv(RESULTS / "x_R1_validation.tsv", sep="\t", index=False)
    print("GATE V:", "PASS" if df["pass"].all() else "FAIL")


def mecfs(n_perm, n_perm_panel):
    rng = np.random.default_rng(SEED)
    cov, spec_cols = load_covar()
    panel = panel_traits()
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t")
    sc_all = {}
    for k_fixed in (None, 3, 6):
        tag = "k_kaiser" if k_fixed is None else f"k{k_fixed}"
        for name in ("decodeme_gwas_1", "decodeme_gwas_2"):
            genes = common_genes(cov, panel + [name])
            Z = np.column_stack([load_genes(x).loc[genes, "ZSTAT"].to_numpy() for x in panel])
            pcs, k, expl = kaiser_pcs(Z, k=k_fixed)
            meta = load_genes(name).loc[genes]
            res, r2 = run_target(meta["ZSTAT"].to_numpy(), meta, pcs, cov, spec_cols, genes, n_perm, rng)
            res["k"] = k
            res["pcs_explain_trait_var"] = r2
            sc_all[(tag, name)] = res
            res.to_csv(RESULTS / f"x_R1_{name}_{tag}.tsv", sep="\t", index=False, float_format="%.5g")
            print(f"{name} {tag}: k={k}, PCs explain {r2:.3f} of ME/CFS gene z; clusters FDR<0.05: {(res['fdr'] < 0.05).sum()}")
            if k_fixed is None and name == "decodeme_gwas_1":
                resid_genes = genes
                Zme = meta["ZSTAT"].to_numpy()
                X = np.column_stack([pcs, np.log(meta["NSNPS"]), np.log(meta["NPARAM"])])
                r, _ = residualise(Zme, X)
                top = pd.DataFrame({"gene": resid_genes, "residual_z": r, "raw_z": Zme}).sort_values("residual_z", ascending=False)
                top.head(200).to_csv(RESULTS / "x_R2_top_residual_genes.tsv", sep="\t", index=False, float_format="%.4g")
    # panel null: each panel trait as target with the other panel traits as PCs
    null_z = {}
    for t in panel:
        others = [x for x in panel if x != t]
        genes = common_genes(cov, others + [t])
        Z = np.column_stack([load_genes(x).loc[genes, "ZSTAT"].to_numpy() for x in others])
        pcs, k, _ = kaiser_pcs(Z)
        meta = load_genes(t).loc[genes]
        res, _ = run_target(meta["ZSTAT"].to_numpy(), meta, pcs, cov, spec_cols, genes, n_perm_panel, rng)
        null_z[t] = calibrate.robust_standardise(res["z"].to_numpy())
        print("panel null done:", t)
    N = np.vstack(list(null_z.values()))
    for name in ("decodeme_gwas_1", "decodeme_gwas_2"):
        res = sc_all[("k_kaiser", name)]
        zs = calibrate.robust_standardise(res["z"].to_numpy())
        res["s"] = (zs - N.mean(0)) / N.std(0)
        res.to_csv(RESULTS / f"x_R1_{name}_k_kaiser.tsv", sep="\t", index=False, float_format="%.5g")
    g1, g2 = sc_all[("k_kaiser", "decodeme_gwas_1")], sc_all[("k_kaiser", "decodeme_gwas_2")]
    call = (g1["fdr"] < 0.05) & (g2["p_emp"] < 0.05) & (g1["s"] >= 2.0)
    out = g1.assign(p_gwas2=g2["p_emp"], called=call)
    out.to_csv(RESULTS / "x_R1_calls.tsv", sep="\t", index=False, float_format="%.5g")
    print("R1 clusters called ME/CFS-specific:", int(call.sum()))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--mecfs", action="store_true")
    ap.add_argument("--perm", type=int, default=2000)
    ap.add_argument("--perm-panel", type=int, default=1000)
    a = ap.parse_args()
    if a.validate:
        validate(a.perm)
    if a.mecfs:
        mecfs(a.perm, a.perm_panel)
