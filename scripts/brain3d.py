"""Map per-cluster enrichment onto brain regions and draw it in 3D (static views and an interactive page).

Region score = cell-weighted mean of the cluster z-values over the cells that came from that region's
dissections (atlas Tissue_* counts). The 3D positions are approximate MNI region centres typed in by hand,
so the figure is a schematic of where the signal sits, not a source-localised map.
"""
import re
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ATLAS = Path("D:/AAN_data/atlas/adult_human_agg.loom")

# approximate MNI centres in mm for the right hemisphere (x > 0); left is the mirror image, midline stays at x = 0
CENTRES = {
    "Amygdala": (24, -4, -20), "Basal forebrain": (8, 6, -10), "Basal nuclei": (20, 6, 4),
    "Hippocampus": (28, -22, -14), "Cerebellum": (24, -62, -32), "Claustrum": (32, 6, 2),
    "Epithalamus": (0, -30, 10), "Extended amygdala": (6, 2, -2), "Hypothalamus": (3, -4, -10),
    "Midbrain": (0, -22, -12), "Medulla": (0, -38, -52), "Paleocortex": (20, 10, -16),
    "Perirhinal cortex": (34, -10, -34), "Pons": (0, -28, -36), "Spinal cord": (0, -44, -70),
    "Thalamus": (10, -18, 8),
}
PREFIX_TO_REGION = [
    (r"^Amygdaloid", "Amygdala"), (r"^Basal forebrain", "Basal forebrain"), (r"^Basal nuclei", "Basal nuclei"),
    (r"(hippocampus|Hippocampus)", "Hippocampus"), (r"^Cerebellum", "Cerebellum"),
    (r"^Cerebral cortex", "Cerebral cortex"), (r"^Claustrum", "Claustrum"), (r"^Epithalamus", "Epithalamus"),
    (r"^Extended amygdala", "Extended amygdala"), (r"^Hypothalamus", "Hypothalamus"),
    (r"^Midbrain", "Midbrain"), (r"^Myelencephalon", "Medulla"), (r"^Paleocortex", "Paleocortex"),
    (r"^Perirhinal", "Perirhinal cortex"), (r"^Pons", "Pons"), (r"^Spinal cord", "Spinal cord"),
    (r"^Thalamus", "Thalamus"),
]


def dissection_region(name):
    for pattern, region in PREFIX_TO_REGION:
        if re.search(pattern, name):
            return region
    return None


def cluster_by_region_counts():
    with h5py.File(ATLAS, "r") as f:
        ids = f["col_attrs"]["Clusters"][:]
        cols = {k[len("Tissue_"):]: f["col_attrs"][k][:] for k in f["col_attrs"].keys() if k.startswith("Tissue_")}
    counts = pd.DataFrame(cols, index=ids)
    regions = {name: dissection_region(name) for name in counts.columns}
    unmapped = [n for n, r in regions.items() if r is None]
    if unmapped:
        print("dissections without a region:", unmapped)
    return counts.T.groupby(pd.Series(regions)).sum().T          # clusters x regions


def region_scores(cluster_z):
    """cluster_z: Series indexed by cluster id."""
    counts = cluster_by_region_counts().loc[cluster_z.index]
    return (counts.mul(cluster_z, axis=0).sum() / counts.sum()).sort_values(ascending=False)


def _points(scores):
    pts = []
    for region, (x, y, z) in CENTRES.items():
        if region not in scores.index:
            continue
        sides = [x] if x == 0 else [x, -x]
        for sx in sides:
            pts.append((region, sx, y, z, scores[region]))
    return pd.DataFrame(pts, columns=["region", "x", "y", "z", "score"])


def _surface():
    from nilearn import datasets, surface
    fs = datasets.fetch_surf_fsaverage("fsaverage5")
    meshes = []
    for hemi in ("left", "right"):
        coords, faces = surface.load_surf_mesh(fs[f"pial_{hemi}"])
        meshes.append((np.asarray(coords), np.asarray(faces)))
    return meshes


def static_views(scores, out_png, title, cmap="YlOrRd"):
    pts = _points(scores)
    vmax = max(scores.max(), 1e-6)
    fig = plt.figure(figsize=(13, 4.6))
    meshes = _surface()
    views = [("lateral", 5, -90), ("dorsal", 90, -90), ("anterior-oblique", 15, 40)]
    for i, (label, elev, azim) in enumerate(views, start=1):
        ax = fig.add_subplot(1, 3, i, projection="3d")
        for coords, faces in meshes:
            ax.plot_trisurf(coords[:, 0], coords[:, 1], faces, coords[:, 2], color="#d9d9d9",
                            alpha=0.10, linewidth=0, shade=False)
        sc = ax.scatter(pts["x"], pts["y"], pts["z"], c=pts["score"], cmap=cmap, vmin=0, vmax=vmax,
                        s=60 + 260 * np.clip(pts["score"], 0, None) / vmax, edgecolor="k", linewidth=0.4,
                        depthshade=False)
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(label, fontsize=9)
        lim = 75
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim - 15, lim); ax.set_zlim(-lim, lim)
    fig.colorbar(sc, ax=fig.axes, shrink=0.6, label="cell-weighted mean enrichment (z)")
    fig.suptitle(title, fontsize=11)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def interactive(scores, out_html, title):
    import plotly.graph_objects as go
    pts = _points(scores)
    fig = go.Figure()
    for coords, faces in _surface():
        fig.add_trace(go.Mesh3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2], i=faces[:, 0], j=faces[:, 1],
                                k=faces[:, 2], color="lightgray", opacity=0.12, hoverinfo="skip", showlegend=False))
    vmax = max(scores.max(), 1e-6)
    fig.add_trace(go.Scatter3d(
        x=pts["x"], y=pts["y"], z=pts["z"], mode="markers", text=pts["region"],
        marker=dict(size=6 + 16 * np.clip(pts["score"], 0, None) / vmax, color=pts["score"],
                    colorscale="YlOrRd", cmin=0, cmax=vmax, colorbar=dict(title="z"), line=dict(width=0.5)),
        hovertemplate="%{text}<br>enrichment z = %{marker.color:.2f}<extra></extra>"))
    fig.update_layout(title=title, scene=dict(aspectmode="data", xaxis_visible=False, yaxis_visible=False,
                                              zaxis_visible=False), margin=dict(l=0, r=0, b=0, t=40))
    fig.write_html(out_html, include_plotlyjs="cdn")
