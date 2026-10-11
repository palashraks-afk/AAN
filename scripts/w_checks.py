"""PREREG v1.14 W1, W3, W4, W5, W6, W7, W8 (read-only on existing results)."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

import brain3d
import calibrate

RES = Path("D:/AAN/results")
DER = Path("D:/AAN_data/derived")
TARGET = "decodeme_gwas_1"
SIX = ["Amex_153", "Amex_175", "DLIT_152", "DLIT_150", "ULIT_121", "Splat_402"]
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


def gsa(path):
    df = pd.read_csv(path, comment="#", sep=r"\s+")
    df = df[df["VARIABLE"].str.match(r"^c\d+$")].copy()
    df["cid"] = df["VARIABLE"].str[1:].astype(int)
    return df.set_index("cid").sort_index()


def main():
    ann = pd.read_csv(DER / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    ids = ann.index[ann["cluster_name"].isin(SIX)]
    nm = ann.loc[ids, "cluster_name"]
    me = calibrate.load_model(TARGET, "A")
    z_me = pd.Series(calibrate.z_from_p(me["P"]), index=me.index)
    zs_me = pd.Series(calibrate.robust_standardise(z_me), index=me.index)
    names = calibrate.panel_names(calibrate.traits_with_results("A"))
    gate = calibrate.passes_heritability_gate(names)
    panel = [n for n in names if gate.get(n) is not False]
    prof = pd.DataFrame({t: calibrate.robust_standardise(calibrate.z_from_p(calibrate.load_model(t, "A")["P"]))
                         for t in panel}, index=me.index)

    def s_of(zs, cols):
        sub = prof[cols]
        return (zs - sub.mean(axis=1)) / sub.std(axis=1, ddof=1)

    # W1
    pt = pd.read_csv(RES / "panel_as_target.tsv", sep="\t").set_index("trait")
    weak = [t for t in panel if t in pt.index and pt.loc[t, "n_sig"] <= 30]
    strong = [t for t in panel if t not in weak]
    say("## W1 matched-power panel")
    say("weak traits (<= 30 sig clusters):", len(weak), weak)
    say("strong traits:", len(strong), strong)
    w1 = pd.DataFrame({"full": s_of(zs_me, panel).loc[ids], "weak_panel": s_of(zs_me, weak).loc[ids],
                       "strong_panel": s_of(zs_me, strong).loc[ids]})
    w1.index = nm.values
    w1["power_robust(s>=2 weak)"] = w1["weak_panel"] >= 2
    say(w1.round(2).to_string())

    # W3
    say("\n## W3 chromosome jackknife of s")
    full_mean, full_sd = prof[panel].mean(axis=1), prof[panel].std(axis=1, ddof=1)
    sj = {}
    for c in range(1, 23):
        p = DER / "magma" / TARGET / f"cluster_loco_chr{c}.gsa.out.txt"
        d = gsa(p)
        zc = calibrate.z_from_p(d["P"])
        zsc = pd.Series(calibrate.robust_standardise(zc), index=d.index)
        sj[c] = ((zsc - full_mean) / full_sd).loc[ids]
    J = pd.DataFrame(sj)
    n = J.shape[1]
    mean = J.mean(axis=1)
    se = np.sqrt((n - 1) / n * ((J.sub(mean, axis=0)) ** 2).sum(axis=1))
    w3 = pd.DataFrame({"s_full": s_of(zs_me, panel).loc[ids], "jack_mean": mean, "jack_se": se,
                       "lo95": mean - 1.96 * se, "hi95": mean + 1.96 * se, "min_loco": J.min(axis=1)})
    w3.index = nm.values
    w3["stable(lo>1.5)"] = w3["lo95"] > 1.5
    say(w3.round(2).to_string())

    # W4
    say("\n## W4 where the cells come from (share of cells by dissection region)")
    cnt = brain3d.cluster_by_region_counts()
    for cid, name in nm.items():
        if name.startswith("Amex"):
            share = (cnt.loc[cid] / cnt.loc[cid].sum()).sort_values(ascending=False)
            say(name, share.head(5).round(3).to_dict())

    # W5
    say("\n## W5 subsets (model A p, Bonferroni 0.05/24 = 0.00208)")
    rows = {}
    for s in ("infectious_onset", "non_infectious_onset", "female", "male"):
        d = pd.read_csv(RES / f"clusters_{TARGET}_{s}.tsv", sep="\t").set_index("cluster_id")
        rows[s] = d.loc[ids, "p_A"]
    w5 = pd.DataFrame(rows)
    w5.index = nm.values
    w5["main"] = me.loc[ids, "P"].values
    say(w5.map(lambda v: f"{v:.2g}").to_string())
    say("passing Bonferroni:", {c: [i for i in w5.index if w5.loc[i, c] < 0.05 / 24] for c in w5.columns if c != "main"})

    # W6
    say("\n## W6 amygdala region score vs panel")
    me_sc = brain3d.region_scores(z_me)
    panel_amy = {}
    for t in panel:
        m = calibrate.load_model(t, "A")
        panel_amy[t] = brain3d.region_scores(pd.Series(calibrate.z_from_p(m["P"]), index=m.index)).get("Amygdala", np.nan)
    pa = pd.Series(panel_amy)
    pct = float((pa < me_sc["Amygdala"]).mean())
    rank = int(me_sc.rank(ascending=False)["Amygdala"])
    say(f"ME/CFS amygdala score {me_sc['Amygdala']:.3f}; rank {rank} of {len(me_sc)} regions; percentile among panel {pct:.2f}")
    say("rule (>90th percentile and top 3):", bool(pct > 0.9 and rank <= 3))
    say("top regions ME/CFS:", me_sc.head(5).round(2).to_dict())

    # W7
    say("\n## W7 power projection (z x sqrt(21620/15579))")
    f = np.sqrt(21620 / 15579)
    for mod in ("A", "B", "C"):
        d = calibrate.load_model(TARGET, mod)
        zp = calibrate.z_from_p(d["P"]) * f
        pp = norm.sf(zp)
        fdr = pd.Series(multipletests(pp, method="fdr_bh")[1], index=d.index)
        base = pd.Series(multipletests(d["P"], method="fdr_bh")[1], index=d.index)
        say(f"model {mod}: clusters FDR<0.05 now {int((base < 0.05).sum())}, projected {int((fdr < 0.05).sum())};",
            "six now/projected:", [(nm[i], round(base[i], 3), round(fdr[i], 3)) for i in ids])

    # W8
    say("\n## W8 other specificity definitions")
    corr = prof[panel].corrwith(zs_me)
    w = (1 - corr).clip(lower=0)
    wm = (prof[panel] * w).sum(axis=1) / w.sum()
    wsd = np.sqrt(((prof[panel].sub(wm, axis=0) ** 2) * w).sum(axis=1) / w.sum())
    s_w = (zs_me - wm) / wsd
    rank_pct = (prof[panel].lt(zs_me, axis=0)).mean(axis=1)
    top_me = zs_me.rank(pct=True) > 0.9
    top_panel = (prof[panel].rank(pct=True) > 0.9).sum(axis=1)
    w8 = pd.DataFrame({"s_weighted": s_w.loc[ids], "rank_pct_vs_panel": rank_pct.loc[ids],
                       "me_in_top10pct": top_me.loc[ids], "panel_traits_top10pct": top_panel.loc[ids]})
    w8.index = nm.values
    w8["(a)rank>=.95"] = w8["rank_pct_vs_panel"] >= 0.95
    w8["(b)s_w>=2"] = w8["s_weighted"] >= 2
    w8["(c)top10 & <=4 of 19 panel"] = w8["me_in_top10pct"] & (w8["panel_traits_top10pct"] <= 4)
    w8["definition_robust"] = w8["(a)rank>=.95"] & w8["(b)s_w>=2"]
    say(w8.round(2).to_string())
    (RES / "w_checks.txt").write_text("\n".join(OUT), encoding="utf-8")


if __name__ == "__main__":
    main()
