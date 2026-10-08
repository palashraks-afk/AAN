"""G3 from the pre-registration: do the positive and negative control traits behave?

    python scripts/check_controls.py
"""
import pandas as pd

import calibrate
from calibrate import EXCITATORY, NON_NEURONAL, load_model, multipletests


def fdr_table(name):
    ann = pd.read_csv(calibrate.DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    m = load_model(name, "A")
    ann = ann.join(m[["P"]])
    ann["fdr"] = multipletests(ann["P"], method="fdr_bh")[1]
    ann["z"] = calibrate.z_from_p(ann["P"])
    return ann


def main():
    checks = []
    scz = fdr_table("schizophrenia")
    ok = ((scz["supercluster"].isin(EXCITATORY)) & (scz["fdr"] < 0.05)).any()
    checks.append(("schizophrenia: an excitatory-neuron cluster at FDR<0.05", ok))

    ad = fdr_table("alzheimer")
    ok = ((ad["supercluster"] == "Microglia") & (ad["fdr"] < 0.05)).any()
    checks.append(("alzheimer: a microglia cluster at FDR<0.05", ok))

    ht = fdr_table("height")
    top3 = ht.sort_values("z", ascending=False).head(3)
    ok = top3["supercluster"].isin(NON_NEURONAL).all()
    checks.append(("height: top-3 clusters all non-neuronal", ok))

    for label, ok in checks:
        print(("PASS  " if ok else "FAIL  ") + label)
    print("G3 overall:", "PASS" if all(ok for _, ok in checks) else "FAIL")
    for name, tab in (("schizophrenia", scz), ("alzheimer", ad), ("height", ht)):
        print("\n", name, "top clusters")
        print(tab.sort_values("z", ascending=False).head(5)[["cluster_name", "supercluster", "z", "fdr"]].to_string())


if __name__ == "__main__":
    main()
