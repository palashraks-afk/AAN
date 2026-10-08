"""3D protein models of candidate genes from the AlphaFold database, drawn as backbone traces coloured by confidence.

    python scripts/structure3d.py ARFGEF2 CSE1L STAU1

For each gene: UniProt accession (reviewed human entry) -> AlphaFold prediction file -> static PNG and interactive HTML.
AlphaFold models are predictions; colour shows the model's own confidence (pLDDT), nothing about disease.
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import requests

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path("D:/AAN/figures/structures")
CACHE = Path("D:/AAN_data/translational/alphafold")


def uniprot_accession(symbol):
    r = requests.get("https://rest.uniprot.org/uniprotkb/search", timeout=60, params={
        "query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true", "fields": "accession", "format": "json",
        "size": 1})
    r.raise_for_status()
    results = r.json().get("results", [])
    return results[0]["primaryAccession"] if results else None


def fetch_model(accession):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"AF-{accession}.pdb"
    if not path.exists():
        meta = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}", timeout=60).json()[0]
        path.write_bytes(requests.get(meta["pdbUrl"], timeout=120).content)
    return path


def read_ca(path):
    xyz, plddt = [], []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            plddt.append(float(line[60:66]))
    return np.array(xyz), np.array(plddt)


def draw(symbol, xyz, plddt):
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(4.6, 4.6))
    ax = fig.add_subplot(111, projection="3d")
    for i in range(len(xyz) - 1):
        ax.plot(xyz[i:i + 2, 0], xyz[i:i + 2, 1], xyz[i:i + 2, 2], color=plt.cm.RdYlBu((plddt[i] - 30) / 70),
                lw=1.1)
    ax.set_axis_off()
    ax.set_title(f"{symbol}  (AlphaFold model, {len(xyz)} residues)", fontsize=9)
    sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=plt.Normalize(30, 100))
    fig.colorbar(sm, ax=ax, shrink=0.55, label="pLDDT confidence")
    fig.savefig(OUT / f"{symbol}.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    import plotly.graph_objects as go
    fig = go.Figure(go.Scatter3d(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], mode="lines+markers",
                                 line=dict(width=4, color=plddt, colorscale="RdYlBu", cmin=30, cmax=100),
                                 marker=dict(size=2, color=plddt, colorscale="RdYlBu", cmin=30, cmax=100),
                                 text=[f"residue {i + 1}, pLDDT {p:.0f}" for i, p in enumerate(plddt)],
                                 hoverinfo="text"))
    fig.update_layout(title=f"{symbol} AlphaFold model", scene=dict(aspectmode="data", xaxis_visible=False,
                                                                    yaxis_visible=False, zaxis_visible=False))
    fig.write_html(OUT / f"{symbol}.html", include_plotlyjs="cdn")


def main(symbols):
    for s in symbols:
        acc = uniprot_accession(s)
        if not acc:
            print(s, "no reviewed UniProt entry")
            continue
        xyz, plddt = read_ca(fetch_model(acc))
        draw(s, xyz, plddt)
        print(s, acc, len(xyz), "residues, mean pLDDT", round(plddt.mean(), 1))


if __name__ == "__main__":
    main(sys.argv[1:])
