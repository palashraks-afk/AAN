"""Inputs for the discovery layer (preregistration/PREREG_v1.2_discovery_layer.md).

    gene_covar_groups.txt   transmitter-class and neuron-subtype groups from the brain atlas   (D1)
    gene_covar_hpa.txt      Human Protein Atlas single-cell types from the whole body          (D2)
    reactome.sets           Reactome pathways as MAGMA gene sets                               (D3)
"""
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

ATLAS = Path("D:/AAN_data/atlas")
DERIVED = Path("D:/AAN_data/derived")
GENELOC = "D:/AAN_data/tools/ref/NCBI38.gene.loc"
MIN_CLUSTERS = 3
MIN_MEAN_CPM = 1.0


def rank_normal(x):
    r = rankdata(x, method="average")
    return norm.ppf((r - 0.5) / len(x))


def transmitter_class(label):
    if pd.isna(label):
        return "unannotated"
    tokens = set(str(label).split())
    glut = any(t.startswith("VGLUT") for t in tokens)
    gaba = "GABA" in tokens
    if gaba and glut:
        return "GABA-glutamate"
    if gaba:
        return "GABAergic"
    if glut and all(t.startswith("VGLUT") for t in tokens):
        return "glutamatergic"
    if "GLY" in tokens:
        return "glycinergic"
    return "monoaminergic-cholinergic-other"


def specificity_table(expr, entrez, names):
    total = expr.sum(axis=1, keepdims=True)
    spec = np.divide(expr, total, out=np.zeros_like(expr, dtype=float), where=total > 0)   # unexpressed genes get 0
    out = pd.DataFrame(np.apply_along_axis(rank_normal, 0, spec), columns=names)
    out.insert(0, "GENE", entrez)
    return out


def atlas_groups():
    covar = pd.read_csv(DERIVED / "gene_covar.txt", sep="\t")       # keeps the same 14,049 genes and avg columns
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    raw = pd.read_excel(ATLAS / "cluster_annotation.xlsx")
    raw.columns = [c.strip() for c in raw.columns]
    raw = raw.set_index("Cluster ID")
    neuron = raw["Class auto-annotation"].eq("NEUR")
    ann["transmitter"] = raw["Neurotransmitter auto-annotation"].map(transmitter_class)
    ann["subtype"] = raw["Subtype auto-annotation"].astype(str).str.split().str[0].replace("nan", np.nan)
    ann = ann[neuron.reindex(ann.index).fillna(False)]

    with h5py.File(ATLAS / "adult_human_agg.loom", "r") as f:
        mat = f["matrix"][:]
        symbols = np.array([g.decode() for g in f["row_attrs"]["Gene"][:]])
        ids = list(f["col_attrs"]["Clusters"][:])
        ncells = pd.Series(f["col_attrs"]["NCells"][:], index=ids)
    cpm = pd.DataFrame(mat / mat.sum(axis=0, keepdims=True) * 1e6, columns=ids)
    loc = pd.read_csv(GENELOC, sep="\t", header=None, names=["entrez", "chr", "s", "e", "st", "symbol"])
    sym2ent = dict(zip(loc["symbol"].drop_duplicates(), loc.drop_duplicates("symbol")["entrez"]))
    cpm["entrez"] = [sym2ent.get(s) for s in symbols]
    cpm = cpm.dropna(subset=["entrez"]).sort_values(ids[0], ascending=False).drop_duplicates("entrez")
    cpm = cpm.set_index(cpm["entrez"].astype(int)).drop(columns="entrez").loc[covar["GENE"].values]

    groups = {}
    for kind in ("transmitter", "subtype"):
        for g, members in ann.groupby(kind).groups.items():
            if kind == "transmitter" and g == "unannotated":
                continue
            if isinstance(g, float) or len(members) < MIN_CLUSTERS:
                continue
            w = ncells.loc[members].to_numpy()
            groups[f"{kind}_{g}"] = (cpm[list(members)].to_numpy() * w).sum(axis=1) / w.sum()
    names = sorted(groups)
    expr = np.column_stack([groups[n] for n in names])
    table = specificity_table(expr, covar["GENE"].values, [n.replace(" ", "_") for n in names])
    table["avg_all"] = covar["avg_all"].values
    table["avg_neuron"] = covar["avg_neuron"].values
    table.to_csv(DERIVED / "gene_covar_groups.txt", sep="\t", index=False, float_format="%.5f")
    ann.assign(n_cells=ncells.loc[ann.index])[["cluster_name", "supercluster", "transmitter", "subtype"]].to_csv(
        DERIVED / "cluster_groups.tsv", sep="\t")
    print("atlas groups:", len(names), names)


def hpa_cell_types():
    z = zipfile.ZipFile(ATLAS / "hpa_rna_single_cell_type.tsv.zip")
    df = pd.read_csv(z.open("rna_single_cell_type.tsv"), sep="\t")
    wide = df.pivot_table(index="Gene name", columns="Cell type", values="nCPM", aggfunc="mean")
    loc = pd.read_csv(GENELOC, sep="\t", header=None, names=["entrez", "chr", "s", "e", "st", "symbol"])
    loc = loc.drop_duplicates("symbol").set_index("symbol")
    wide = wide[wide.index.isin(loc.index)]
    wide = wide[wide.mean(axis=1) >= MIN_MEAN_CPM]
    entrez = loc.loc[wide.index, "entrez"].to_numpy()
    table = specificity_table(wide.to_numpy(), entrez, [c.replace(" ", "_") for c in wide.columns])
    table["avg_all"] = np.log10(wide.mean(axis=1).to_numpy() + 1)
    table.to_csv(DERIVED / "gene_covar_hpa.txt", sep="\t", index=False, float_format="%.5f")
    print("HPA cell types:", wide.shape[1], "genes:", wide.shape[0])


def reactome_sets():
    z = zipfile.ZipFile("D:/AAN_data/translational/ReactomePathways.gmt.zip")
    loc = pd.read_csv(GENELOC, sep="\t", header=None, names=["entrez", "chr", "s", "e", "st", "symbol"])
    sym2ent = dict(zip(loc["symbol"], loc["entrez"]))
    kept, names = 0, []
    with open(DERIVED / "reactome.sets", "w") as out:
        for line in z.read("ReactomePathways.gmt").decode().splitlines():
            parts = line.split("\t")
            title, rid, genes = parts[0], parts[1], parts[2:]
            ents = sorted({str(sym2ent[g]) for g in genes if g in sym2ent})
            if 10 <= len(ents) <= 500:
                out.write(rid + " " + " ".join(ents) + "\n")
                names.append((rid, title, len(ents)))
                kept += 1
    pd.DataFrame(names, columns=["set", "title", "n_genes"]).to_csv(DERIVED / "reactome_names.tsv", sep="\t", index=False)
    print("Reactome sets kept:", kept)


if __name__ == "__main__":
    atlas_groups()
    hpa_cell_types()
    reactome_sets()
