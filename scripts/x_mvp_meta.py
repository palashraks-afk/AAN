"""PREREG v1.15 M1 + v1.16: MVP-only ME/CFS GWAS (independent) and the DecodeME+MVP meta-analysis (descriptive only).

    python scripts/x_mvp_meta.py
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import calibrate

HERE = Path(__file__).resolve().parent
DERIVED = Path("D:/AAN_data/derived")
RESULTS = Path("D:/AAN/results")
MVP = Path("D:/AAN_data/mvp/GCST90479178.tsv.gz")
META = Path("D:/AAN_data/maccallini/GWAS_METAL_DME_1_MVP_GRCh38.tsv.gz")
META_H = Path("D:/AAN_data/maccallini/meta_harmonised.tsv.gz")
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]
ELEVEN = SIX + ["DLIT_136", "DLIT_151", "ULIT_126", "DLIT_148", "DLIT_147"]


def run(cmd):
    print(" ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], cwd=HERE).returncode == 0


def convert_meta():
    if META_H.exists():
        return
    m = pd.read_csv(META, sep="\t", usecols=["SNP", "CHR", "BP", "N", "P"], dtype={"CHR": str})
    out = pd.DataFrame({"chromosome": m["CHR"], "base_pair_location": m["BP"], "rsid": m["SNP"], "p_value": m["P"], "n": m["N"]})
    out.to_csv(META_H, sep="\t", index=False, compression="gzip")


def analyse(name, path, n_total):
    f = DERIVED / "ldsc_results.tsv"
    have = f.exists() and name in pd.read_csv(f, sep="\t")["name"].values
    if not have:
        run([sys.executable, "-u", "ldsc.py", "--name", name, "--format", "harmonised", "--gwas", path, "--n", n_total])
    t = pd.read_csv(f, sep="\t").drop_duplicates("name", keep="last").set_index("name")
    z_h2 = float(t.loc[name, "h2_z"])
    if not any((DERIVED / "magma" / name / f"cluster_A.gsa.out{s}").exists() for s in (".txt", "")):
        run([sys.executable, "-u", "run_magma.py", "--name", name, "--format", "harmonised", "--gwas", path, "--n", n_total, "--jobs", 8])
    m = calibrate.load_model(name, "A")
    m["z"] = calibrate.z_from_p(m["P"].to_numpy())
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    m = m.join(ann[["cluster_name"]])
    six = m[m["cluster_name"].isin(SIX)]
    obs = six["z"].mean()
    rng = np.random.default_rng(20261010)
    null = np.array([rng.choice(m["z"].to_numpy(), 6, replace=False).mean() for _ in range(10000)])
    p = (1 + (null >= obs).sum()) / 10001
    return m, {"trait": name, "h2_z": z_h2, "mean_z_six": obs, "null_95th": np.percentile(null, 95), "p": p, "n_positive_of_6": int((six["z"] > 0).sum())}


def main():
    convert_meta()
    rows, tabs = [], {}
    for name, path, n in (("mvp_mecfs", MVP, 443093), ("decodeme_mvp_meta", META_H, 700000)):
        m, r = analyse(name, path, n)
        rows.append(r)
        tabs[name] = m.set_index("cluster_name")["z"]
    dm = calibrate.load_model("decodeme_gwas_1", "A")
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    dm["z"] = calibrate.z_from_p(dm["P"].to_numpy())
    tabs["decodeme_gwas_1"] = dm.join(ann[["cluster_name"]]).set_index("cluster_name")["z"]
    six = pd.DataFrame({k: v.reindex(ELEVEN) for k, v in tabs.items()}).round(2)
    six.to_csv(RESULTS / "x_M1_cluster_z_by_cohort.tsv", sep="\t")
    res = pd.DataFrame(rows)
    res["call"] = np.where(res["h2_z"] < 4, "below heritability gate: descriptive only", np.where(res["p"] < 0.05, "supported", "not supported"))
    res.loc[res["trait"] == "decodeme_mvp_meta", "call"] = "descriptive only (contains DecodeME cases)"
    res.to_csv(RESULTS / "x_M1_mvp_meta.tsv", sep="\t", index=False, float_format="%.4g")
    print(res.to_string(index=False))
    print(six.to_string())


if __name__ == "__main__":
    main()
