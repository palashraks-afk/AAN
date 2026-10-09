"""L3 (PREREG v1.5): a 3D map of the 461 brain cell clusters, coloured by a trait's enrichment.

Clusters are placed by the first three principal components of their log expression, so clusters of the same family sit
together. This is a picture of cell similarity, not of anatomy; the colour is the only thing that comes from the GWAS.

    python scripts/landscape3d.py --trait decodeme_gwas_1
"""
import argparse
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pandas as pd

import calibrate

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ATLAS = Path("D:/AAN_data/atlas/adult_human_agg.loom")
DERIVED = Path("D:/AAN_data/derived")
FIG = Path("D:/AAN/figures")


def cluster_coordinates():
    with h5py.File(ATLAS, "r") as f:
        mat = f["matrix"][:]
        ids = list(f["col_attrs"]["Clusters"][:])
        ncells = pd.Series(f["col_attrs"]["NCells"][:], index=ids)
    cpm = mat / mat.sum(axis=0, keepdims=True) * 1e6
    keep = cpm.mean(axis=1) >= 1.0
    x = np.log10(cpm[keep] + 1)
    x = x - x.mean(axis=1, keepdims=True)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    coords = (vt[:3].T * s[:3])
    var = (s[:3] ** 2) / (s ** 2).sum()
    return pd.DataFrame(coords, index=ids, columns=["pc1", "pc2", "pc3"]), ncells, var


def draw(trait, values, label, tag):
    coords, ncells, var = cluster_coordinates()
    ann = pd.read_csv(DERIVED / "cluster_annotation.tsv", sep="\t").set_index("cluster_id")
    df = coords.join(ann[["cluster_name", "supercluster"]]).join(ncells.rename("n_cells"))
    df["value"] = values.reindex(df.index)
    vmax = float(np.nanpercentile(df["value"], 99))

    fig = plt.figure(figsize=(13, 4.6))
    for i, (elev, azim) in enumerate([(20, 30), (20, 120), (80, 270)], start=1):
        ax = fig.add_subplot(1, 3, i, projection="3d")
        sc = ax.scatter(df["pc1"], df["pc2"], df["pc3"], c=df["value"], cmap="YlOrRd", vmin=0, vmax=vmax,
                        s=6 + 40 * np.log10(df["n_cells"].clip(lower=10)) / 4, edgecolor="none", depthshade=False)
        ax.view_init(elev=elev, azim=azim)
        ax.set_xlabel(f"PC1 ({var[0]:.0%})", fontsize=7); ax.set_ylabel(f"PC2 ({var[1]:.0%})", fontsize=7)
        ax.set_zlabel(f"PC3 ({var[2]:.0%})", fontsize=7)
        ax.tick_params(labelsize=5)
    fig.colorbar(sc, ax=fig.axes, shrink=0.6, label=label)
    fig.suptitle(f"{trait}: enrichment across the brain cell-type landscape", fontsize=10)
    fig.savefig(FIG / f"fig8_landscape_{tag}.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    import plotly.graph_objects as go
    top = df.nlargest(8, "value")
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=df["pc1"], y=df["pc2"], z=df["pc3"], mode="markers",
        marker=dict(size=3 + 8 * np.log10(df["n_cells"].clip(lower=10)) / 4, color=df["value"], colorscale="YlOrRd",
                    cmin=0, cmax=vmax, colorbar=dict(title=label), line=dict(width=0)),
        text=[f"{n} ({s})<br>{label} = {v:.2f}" for n, s, v in zip(df["cluster_name"], df["supercluster"], df["value"])],
        hoverinfo="text"))
    fig.add_trace(go.Scatter3d(x=top["pc1"], y=top["pc2"], z=top["pc3"], mode="text", text=top["cluster_name"],
                               textfont=dict(size=9), showlegend=False, hoverinfo="skip"))
    fig.update_layout(title=f"{trait}: brain cell-type landscape", margin=dict(l=0, r=0, b=0, t=40),
                      scene=dict(xaxis_title=f"PC1 ({var[0]:.0%})", yaxis_title=f"PC2 ({var[1]:.0%})",
                                 zaxis_title=f"PC3 ({var[2]:.0%})"))
    fig.write_html(FIG / f"fig8_landscape_{tag}.html", include_plotlyjs="cdn")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True)
    a = ap.parse_args()
    FIG.mkdir(exist_ok=True)
    model = calibrate.load_model(a.trait, "A")
    z = pd.Series(calibrate.z_from_p(model["P"].to_numpy()), index=model.index)
    draw(a.trait, z, "enrichment z", a.trait)
