"""PREREG v1.12: use each panel trait as the target, with the other panel traits as its panel."""
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import calibrate

RES = Path("D:/AAN/results")
EXPECT = {
    "alzheimer": ["Microglia"],
    "rheumatoid_arthritis": ["Microglia", "Miscellaneous"],
    "ibd": ["Microglia", "Miscellaneous"],
    "schizophrenia": [],      # any excitatory-neuron supercluster, matched by name below
    "height": ["Vascular", "Miscellaneous"],
}


def main():
    ann = pd.read_csv("D:/AAN_data/derived/cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    names = calibrate.panel_names(calibrate.traits_with_results("A"))
    gate = calibrate.passes_heritability_gate(names)
    panel = [n for n in names if gate.get(n) is not False]
    models = {t: calibrate.load_model(t, "A") for t in panel}
    z = pd.DataFrame({t: calibrate.z_from_p(models[t]["P"]) for t in panel}, index=ann.index)
    prof = z.apply(calibrate.robust_standardise)
    rows, spec_super = [], {}
    for t in panel:
        others = [p for p in panel if p != t]
        s = (prof[t] - prof[others].mean(axis=1)) / prof[others].std(axis=1, ddof=1)
        fdr = pd.Series(multipletests(models[t]["P"], method="fdr_bh")[1], index=ann.index)
        sig, spec = fdr < 0.05, (fdr < 0.05) & (s >= 2)
        sup = ann.loc[spec, "supercluster"]
        for u in sup.unique():
            spec_super.setdefault(u, set()).add(t)
        exp = EXPECT.get(t)
        if t == "schizophrenia":
            ok = bool(sup.str.contains("excitatory|intratelencephalic|Deep-layer|Upper-layer", case=False).any())
        elif exp:
            ok = bool(sup.isin(exp).any())
        else:
            ok = None
        rows.append({"trait": t, "n_sig": int(sig.sum()), "n_specific": int(spec.sum()),
                     "share_specific": spec.sum() / sig.sum() if sig.sum() else np.nan,
                     "expected_recovered": ok, "specific_superclusters": ";".join(sorted(sup.unique()))[:120]})
    out = pd.DataFrame(rows)
    out.to_csv(RES / "panel_as_target.tsv", sep="\t", index=False, float_format="%.3f")
    print(out.to_string(index=False))
    have = out[out["expected_recovered"].notna()]
    print("\nexpected-biology recovered:", int(have["expected_recovered"].sum()), "of", len(have))
    sh = out.loc[out.n_sig >= 3, "share_specific"]
    print("traits with >= 3 significant clusters:", len(sh), "| median share specific:", round(sh.median(), 3),
          "| 90th percentile:", round(sh.quantile(0.9), 3), "| ME/CFS share: 6/11 = 0.545")
    cls = pd.Series({k: len(v) for k, v in spec_super.items()}).sort_values(ascending=False)
    print("\nclasses specific in most traits:\n", cls.head(8).to_string())


if __name__ == "__main__":
    main()
