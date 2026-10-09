"""Run the discovery-layer MAGMA tests for every trait that already has merged gene results.

Safe to restart: skips outputs that exist. Reactome is run for the ME/CFS GWAS only (plus the top-cluster
conditioning for gwas_1 and gwas_2, using the cluster with the smallest model-A p, see PREREG v1.2).
"""
from pathlib import Path

import calibrate
import gene_property

DERIVED = Path("D:/AAN_data/derived")


def have(name, prefix):
    base = DERIVED / "magma" / name
    return any((base / f"{prefix}.gsa.out{s}").exists() for s in (".txt", ""))


def main():
    for d in sorted((DERIVED / "magma").iterdir()):
        name = d.name
        if not (d / "genes.genes.raw").exists():
            continue
        if not have(name, "groups_B"):
            gene_property.run_groups(name)
        if not have(name, "hpa_A"):
            gene_property.run_hpa(name)
        if not have(name, "regions_A") and (DERIVED / "gene_covar_regions.txt").exists():
            gene_property.run_regions(name)
        if name.startswith("decodeme_gwas") and name in ("decodeme_gwas_1", "decodeme_gwas_2"):
            if not have(name, "reactome_A"):
                gene_property.run_reactome(name)
            if not have(name, "reactome_C") and (d / "cluster_A.gsa.out.txt").exists() | (d / "cluster_A.gsa.out").exists():
                top = int(calibrate.load_model("decodeme_gwas_1", "A")["P"].idxmin())
                print("top cluster for pathway conditioning:", top)
                gene_property.run_reactome(name, top)


if __name__ == "__main__":
    main()
