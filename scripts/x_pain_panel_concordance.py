"""PREREG v1.18: pain-augmented panel sensitivity (P1) and Mac/Windows concordance (P2).

    python scripts/x_pain_panel_concordance.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import calibrate
import x_specific_component as x

DERIVED = Path("D:/AAN_data/derived")
RESULTS = Path("D:/AAN/results")
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]


def zprof(name):
    m = calibrate.load_model(name, "A")
    return m["P"], calibrate.z_from_p(m["P"].to_numpy())


def s_scores(target, panel):
    _, zt = zprof(target)
    zt = calibrate.robust_standardise(zt)
    N = np.vstack([calibrate.robust_standardise(zprof(t)[1]) for t in panel])
    return (zt - N.mean(0)) / N.std(0)


def main():
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    panel = x.panel_traits()
    s19 = s_scores("decodeme_gwas_1", panel)
    s21 = s_scores("decodeme_gwas_1", panel + ["finngen_pain", "finngen_fibromyalgia"])
    s20p = s_scores("decodeme_gwas_1", panel + ["finngen_pain"])
    names = ann["cluster_name"].to_numpy()
    out = pd.DataFrame({"cluster": names, "s_panel19": s19, "s_plus_pain": s20p, "s_plus_pain_fibro": s21}, index=ann.index)
    six = out[out["cluster"].isin(SIX)].round(2)
    n_stay = int((six["s_plus_pain_fibro"] >= 2).sum())
    out.to_csv(RESULTS / "x_P1_pain_augmented_panel.tsv", sep="\t", float_format="%.4g")
    print("P1 specificity score of the six clusters (panel of 19 vs plus pain vs plus pain and fibromyalgia):")
    print(six.to_string(index=False))
    print("stay specific (s >= 2) in the 21-trait panel:", n_stay, "of 6")
    # direct comparison: each cluster's robustly standardised z in ME/CFS versus pain, fibromyalgia
    zs = {t: calibrate.robust_standardise(zprof(t)[1]) for t in ("decodeme_gwas_1", "finngen_pain", "finngen_fibromyalgia")}
    cmp_ = pd.DataFrame({"cluster": names, **{k: v for k, v in zs.items()}}, index=ann.index)
    six_cmp = cmp_[cmp_["cluster"].isin(SIX)].round(2)
    print("robust z within trait (ME/CFS vs pain vs fibromyalgia) for the six:")
    print(six_cmp.to_string(index=False))
    print("Spearman across all 461 clusters: ME/CFS vs pain =", round(spearmanr(zs["decodeme_gwas_1"], zs["finngen_pain"])[0], 3), "; ME/CFS vs fibromyalgia =", round(spearmanr(zs["decodeme_gwas_1"], zs["finngen_fibromyalgia"])[0], 3))
    cmp_.to_csv(RESULTS / "x_P1_mecfs_vs_pain_z.tsv", sep="\t", float_format="%.4g")
    # P2: concordance with the Windows pipeline
    win = pd.read_csv(RESULTS / "clusters_decodeme_gwas_1.tsv", sep="\t").set_index("cluster_id")
    p_mac, z_mac = zprof("decodeme_gwas_1")
    lp_mac = -np.log10(p_mac.reindex(win.index))
    lp_win = -np.log10(win["p_A"])
    r_p = spearmanr(lp_mac, lp_win)[0]
    r_s = spearmanr(out.loc[win.index, "s_panel19"], win["s"])[0]
    r_z = np.corrcoef(lp_mac, lp_win)[0, 1]
    print(f"P2 Mac vs Windows over 461 clusters: Spearman(-log10 p_A) = {r_p:.4f}, Pearson = {r_z:.4f}, Spearman(s) = {r_s:.4f}")
    top = win.sort_values("p_A").head(11).index
    print("the 11 FDR-significant clusters on the Mac: p_A <", float(p_mac.reindex(top).max()), "; Windows max p_A", float(win.loc[top, "p_A"].max()))
    pd.DataFrame([{"spearman_log10p": r_p, "pearson_log10p": r_z, "spearman_s": r_s, "six_stay_specific_in_21_trait_panel": n_stay}]).to_csv(RESULTS / "x_P2_concordance.tsv", sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
