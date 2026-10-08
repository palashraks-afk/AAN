"""Readers that turn different GWAS file layouts into one table: rsid, chrom, pos, p, n.

Positions stay in whatever build the file uses; the caller picks the matching gene file.
"""
import numpy as np
import pandas as pd

AUTOSOMES = {str(i) for i in range(1, 23)}


def read_decodeme(path, rsid_map, min_maf=0.01):
    cols = ["CHROM", "GENPOS", "ID", "A1FREQ", "N", "LOG10P"]
    g = pd.read_csv(path, sep=r"\s+", usecols=cols, dtype={"CHROM": str})
    maf = np.minimum(g["A1FREQ"], 1 - g["A1FREQ"])
    g = g[maf >= min_maf]
    g = g.merge(rsid_map[["ID", "rsid"]], on="ID", how="inner")
    g["p"] = 10.0 ** (-g["LOG10P"])
    out = g.rename(columns={"CHROM": "chrom", "GENPOS": "pos", "N": "n"})
    return out[["rsid", "chrom", "pos", "p", "n"]]


def read_harmonised(path, total_n):
    """EBI GWAS Catalog harmonised file (GRCh38). Uses hm_rsid / hm_chrom / hm_pos."""
    head = pd.read_csv(path, sep="\t", nrows=1)
    want = ["hm_rsid", "hm_chrom", "hm_pos", "p_value"]
    extra = [c for c in ("n", "n_cas", "n_con") if c in head.columns]
    g = pd.read_csv(path, sep="\t", usecols=want + extra, dtype={"hm_chrom": str})
    g = g.dropna(subset=["hm_rsid", "hm_pos", "p_value"])
    g = g[g["hm_chrom"].isin(AUTOSOMES)]
    g = g[g["p_value"] > 0]
    if "n" in extra and g["n"].notna().any():
        n = g["n"].fillna(total_n)
    elif {"n_cas", "n_con"} <= set(extra) and g["n_cas"].notna().any():
        n = (g["n_cas"] + g["n_con"]).fillna(total_n)
    else:
        n = total_n
    out = pd.DataFrame({"rsid": g["hm_rsid"], "chrom": g["hm_chrom"], "pos": g["hm_pos"].astype(int),
                        "p": g["p_value"], "n": n})
    return out.drop_duplicates("rsid")


def read_pgc3(path):
    """PGC3 schizophrenia release (GRCh37), columns CHROM ID POS ... PVAL NCAS NCON NEFF."""
    skip = 0
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            skip += 1
    g = pd.read_csv(path, sep="\t", skiprows=skip, dtype={"CHROM": str}, encoding_errors="replace",
                    usecols=["CHROM", "ID", "POS", "PVAL", "NCAS", "NCON"])
    g = g[g["CHROM"].isin(AUTOSOMES) & g["ID"].astype(str).str.startswith("rs")]
    out = pd.DataFrame({"rsid": g["ID"], "chrom": g["CHROM"], "pos": g["POS"],
                        "p": g["PVAL"], "n": g["NCAS"] + g["NCON"]})
    return out[out["p"] > 0].drop_duplicates("rsid")


def open_text(path):
    import gzip
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")
