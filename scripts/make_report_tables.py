"""Markdown tables for paper/REPORT.md, built from results/ so the numbers cannot be mistyped."""
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

import calibrate

RES = Path("D:/AAN/results")
DER = Path("D:/AAN_data/derived")
OUT = Path("D:/AAN/paper/tables.md")
ORDER_SIX = ["Amex_153", "Amex_175", "ULIT_121", "Splat_402", "DLIT_152", "DLIT_150"]


def md(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(v) for v in r.values) + " |")
    return "\n".join(lines)


def fmt_p(p):
    return f"{p:.1e}" if p < 0.001 else f"{p:.3f}"


def main():
    c1 = pd.read_csv(RES / "clusters_decodeme_gwas_1.tsv", sep="\t").set_index("cluster_id")
    c2 = pd.read_csv(RES / "clusters_decodeme_gwas_2.tsv", sep="\t").set_index("cluster_id")
    sl = pd.read_csv(RES / "sldsc_decodeme_gwas_1.tsv", sep="\t")
    sl["cid"] = sl["cluster"].str[1:].astype(int)
    sl = sl.set_index("cid")
    pm = pd.read_csv(RES / "permutation_decodeme_gwas_1.tsv", sep="\t").set_index("cluster")
    mc = calibrate.load_model("decodeme_gwas_1", "C")
    mc_fdr = pd.Series(multipletests(mc["P"], method="fdr_bh")[1], index=mc.index)

    sig = c1[c1["fdr_A"] < 0.05].sort_values("p_A")
    rows = []
    for cid, r in sig.iterrows():
        rows.append({"Cluster": r["cluster_name"], "Cell class": r["supercluster"], "z": f"{r['z']:.2f}", "FDR (model A)": f"{r['fdr_A']:.3f}",
                     "p (model B)": fmt_p(r["p_B"]), "FDR (model C)": f"{mc_fdr[cid]:.3f}", "s": f"{r['s']:.2f}",
                     "Specific": "yes" if r["T2_specific"] else "no",
                     "S-LDSC p": fmt_p(sl.loc[cid, "p_one_sided"]), "Permutation p": fmt_p(pm.loc[cid, "empirical_p"]),
                     "2nd analysis p": fmt_p(c2.loc[cid, "p_A"])})
    t1 = md(pd.DataFrame(rows))

    pr = pd.read_csv(RES / "panel_robustness.tsv", sep="\t").set_index("cluster").loc[ORDER_SIX]
    ls = {}
    for line in (RES / "w_checks.txt").read_text(encoding="utf-8").splitlines():
        pass
    t3 = md(pd.DataFrame({"Cluster": pr.index, "s (all 19)": pr["full panel"].round(2).values, "Lowest s, drop one trait": pr["R2 min s"].round(2).values,
                          "Share of half panels with s >= 2": (pr["R3 share >= 2"] * 100).round(0).astype(int).astype(str).add("%").values,
                          "s without similar traits": pr["R4 s without similar"].round(2).values,
                          "s without psychiatric traits": pr["R4 s without psych"].round(2).values}))

    pt = pd.read_csv(RES / "panel_as_target.tsv", sep="\t")
    pt = pt[pt["n_sig"] >= 3].sort_values("n_sig", ascending=False)
    t4 = md(pd.DataFrame({"Trait": pt["trait"], "Enriched clusters": pt["n_sig"], "Specific": pt["n_specific"],
                          "Share specific": (pt["share_specific"] * 100).round(1).astype(str).add("%"),
                          "Specific cell classes": pt["specific_superclusters"].fillna("none").str.slice(0, 60)}))

    panel = pd.read_csv("D:/AAN/data/panel.tsv", sep="\t")
    keep = set(calibrate.panel_names(calibrate.traits_with_results("A")))
    gate = calibrate.passes_heritability_gate(sorted(keep))
    used = [n for n in keep if gate.get(n) is not False]
    panel = panel[panel["label"].isin(used)].sort_values("label")
    t5 = md(pd.DataFrame({"Trait": panel["trait"], "Study": panel["first_author"] + " " + panel["year"].astype(str),
                          "Sample size": panel["n_total"].map(lambda v: f"{int(v):,}"), "Accession": panel["accession"]}))

    cg = pd.read_csv(RES / "candidate_genes_decodeme_gwas_1.tsv", sep="\t").head(10)
    t6 = md(pd.DataFrame({"Gene": cg["symbol"], "Gene z": cg["ZSTAT"].round(2), "Gene p": cg["P"].map(fmt_p), "Rank percentile": cg["rank_pct"].round(3)}))

    OUT.write_text("\n\n".join([f"<!--T1-->\n{t1}", f"<!--T3-->\n{t3}", f"<!--T4-->\n{t4}", f"<!--T5-->\n{t5}", f"<!--T6-->\n{t6}"]), encoding="utf-8")
    print(OUT.read_text(encoding="utf-8")[:3500])


if __name__ == "__main__":
    main()
