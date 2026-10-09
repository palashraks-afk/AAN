"""Readers that turn different GWAS file layouts into one table: rsid, chrom, pos, p, n.

Positions stay in whatever build the file uses; the caller picks the matching gene file.
"""
import numpy as np
import pandas as pd

AUTOSOMES = {str(i) for i in range(1, 23)}


QC_LIST = "D:/AAN_data/decodeme/gwas_qced.var.gz"


def read_decodeme(path, rsid_map, min_maf=0.01, qc_list=QC_LIST):
    """DecodeME REGENIE output, keeping only the variants on the team's QC list.

    The shared files are raw REGENIE output (filtered on MAF only). The team's README says to keep
    the variants that passed QC (INFO >= 0.9 plus a DENTIST-style test); without that filter
    there are hundreds of artefact SNPs at p < 1e-30 in the pericentromeric regions of chr15 and chr21.
    """
    cols = ["CHROM", "GENPOS", "ID", "A1FREQ", "N", "LOG10P"]
    g = pd.read_csv(path, sep=r"\s+", usecols=cols, dtype={"CHROM": str})
    if qc_list:
        passed = pd.read_csv(qc_list, header=None, names=["ID"])["ID"]
        g = g[g["ID"].isin(set(passed))]
    maf = np.minimum(g["A1FREQ"], 1 - g["A1FREQ"])
    g = g[maf >= min_maf]
    g = g.merge(rsid_map[["ID", "rsid"]], on="ID", how="inner")
    g["p"] = 10.0 ** (-g["LOG10P"])
    out = g.rename(columns={"CHROM": "chrom", "GENPOS": "pos", "N": "n"})
    return out[["rsid", "chrom", "pos", "p", "n"]]


def read_harmonised(path, total_n):
    """EBI GWAS Catalog harmonised file (GRCh38). Uses hm_rsid / hm_chrom / hm_pos."""
    # two layouts exist among the harmonised files: hm_* columns, or the newer chromosome/base_pair_location/rsid layout
    # (both are GRCh38 after harmonisation)
    head = pd.read_csv(path, sep="\t", nrows=1)
    cols = set(head.columns)
    chrom_c = "hm_chrom" if "hm_chrom" in cols else "chromosome"
    pos_c = "hm_pos" if "hm_pos" in cols else "base_pair_location"
    rs_c = "hm_rsid" if "hm_rsid" in cols else ("rsid" if "rsid" in cols else "variant_id")
    n_cols = [c for c in ("n", "N", "n_cas", "n_con", "N_cases", "N_controls") if c in cols]
    g = pd.read_csv(path, sep="\t", usecols=[chrom_c, pos_c, rs_c, "p_value"] + n_cols, dtype={chrom_c: str})
    g = g.dropna(subset=[rs_c, pos_c, "p_value"])
    g = g[g[chrom_c].isin(AUTOSOMES) & (g["p_value"] > 0)]
    g = g[g[rs_c].astype(str).str.startswith("rs")]
    if "n" in n_cols and g["n"].notna().any():
        n = g["n"].fillna(total_n)
    elif {"n_cas", "n_con"} <= set(n_cols) and g["n_cas"].notna().any():
        n = (g["n_cas"] + g["n_con"]).fillna(total_n)
    elif {"N_cases", "N_controls"} <= set(n_cols) and g["N_cases"].notna().any():
        n = (g["N_cases"] + g["N_controls"]).fillna(total_n)
    else:
        n = total_n
    out = pd.DataFrame({"rsid": g[rs_c], "chrom": g[chrom_c], "pos": g[pos_c].astype(int), "p": g["p_value"], "n": n})
    out = out.drop_duplicates("rsid")
    if len(out) < 300_000:
        raise ValueError(f"only {len(out)} usable SNPs in {path}; check the file layout")
    return out


def read_finngen(path, total_n, min_maf=0.01):
    """FinnGen R13 summary statistics (GRCh38, rsIDs included, no per-SNP N)."""
    cols = ["#chrom", "pos", "rsids", "pval", "af_alt"]
    g = pd.read_csv(path, sep="\t", usecols=cols, dtype={"#chrom": str})
    g = g[g["#chrom"].isin(AUTOSOMES) & g["rsids"].notna() & (g["pval"] > 0)]
    maf = np.minimum(g["af_alt"], 1 - g["af_alt"])
    g = g[maf >= min_maf]
    g["rsid"] = g["rsids"].str.split(",").str[0]
    g = g[g["rsid"].str.startswith("rs")]
    out = pd.DataFrame({"rsid": g["rsid"], "chrom": g["#chrom"], "pos": g["pos"].astype(int),
                        "p": g["pval"], "n": total_n})
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
