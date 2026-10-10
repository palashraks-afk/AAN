"""R2 to R4 of preregistration/PREREG_v1.8_robustness.md: how much does the specificity score s depend on the choice of panel?

    python scripts/panel_robustness.py

Recomputes s for the clusters declared specific in DISCOVERY.md with panels that leave traits out. Reads only files that already exist.
"""
from pathlib import Path

import numpy as np
import pandas as pd

import calibrate

RES = Path("D:/AAN/results")
TARGET = "decodeme_gwas_1"
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]
CORRELATED = ["schizophrenia", "education_years", "bmi", "neuroticism"]
PSYCH = ["depression_broad", "anxiety", "neuroticism", "adhd", "autism", "education_years"]
N_HALF = 1000
SEED = 20261010


def s_scores(z_std, profiles, panel):
    sub = profiles[panel]
    return (z_std - sub.mean(axis=1)) / sub.std(axis=1, ddof=1)


def main():
    ann = pd.read_csv("D:/AAN_data/derived/cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    me = calibrate.load_model(TARGET, "A")
    z_std = pd.Series(calibrate.robust_standardise(calibrate.z_from_p(me["P"])), index=me.index)

    names = calibrate.panel_names(calibrate.traits_with_results("A"))
    gate = calibrate.passes_heritability_gate(names)
    panel = [n for n in names if gate.get(n) is not False]
    profiles = pd.DataFrame({t: calibrate.robust_standardise(calibrate.z_from_p(calibrate.load_model(t, "A")["P"]))
                             for t in panel}, index=me.index)
    ids = ann.index[ann["cluster_name"].isin(SIX)]
    names_of = ann.loc[ids, "cluster_name"]
    print("panel traits:", len(panel))

    full = s_scores(z_std, profiles, panel).loc[ids]
    rows = {"full panel": full}

    # R2 leave one trait out
    loo = pd.DataFrame({t: s_scores(z_std, profiles, [p for p in panel if p != t]).loc[ids] for t in panel})
    rows["R2 min s"] = loo.min(axis=1)
    rows["R2 share >= 2"] = (loo >= 2).mean(axis=1)
    worst = loo.idxmin(axis=1)

    # R3 random half panels
    rng = np.random.default_rng(SEED)
    k = len(panel) // 2
    hits = pd.Series(0.0, index=ids)
    for _ in range(N_HALF):
        pick = list(rng.choice(panel, size=k, replace=False))
        hits += (s_scores(z_std, profiles, pick).loc[ids] >= 2).astype(float)
    rows["R3 share >= 2"] = hits / N_HALF

    # R4 without the correlated or psychiatric traits
    rows["R4 s without similar"] = s_scores(z_std, profiles, [p for p in panel if p not in CORRELATED]).loc[ids]
    rows["R4 s without psych"] = s_scores(z_std, profiles, [p for p in panel if p not in PSYCH]).loc[ids]

    out = pd.DataFrame(rows)
    out.insert(0, "cluster", names_of)
    out["R2 worst trait left out"] = worst
    out["robust"] = (out["R2 share >= 2"] >= 0.9) & (out["R3 share >= 2"] >= 0.8)
    out.to_csv(RES / "panel_robustness.tsv", sep="\t", index=False, float_format="%.3f")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
