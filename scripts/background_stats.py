"""Background numbers for the report: how much research ME/CFS gets compared with similar-sized neurological conditions.

Publication counts come from PubMed (NCBI E-utilities), so they can be rerun. Patient numbers and funding are taken
from the sources listed in docs/08_background_and_prior_models.md and are typed in below with their citations.

    python scripts/background_stats.py
"""
import time
from pathlib import Path

import matplotlib
import pandas as pd
import requests

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
OUT = Path("D:/AAN")
YEARS = list(range(2005, 2026))

TERMS = {
    "ME/CFS": '("myalgic encephalomyelitis"[Title/Abstract] OR "chronic fatigue syndrome"[Title/Abstract])',
    "Multiple sclerosis": '"multiple sclerosis"[Title/Abstract]',
    "Parkinson's disease": '"Parkinson disease"[Title/Abstract] OR "Parkinson\'s disease"[Title/Abstract]',
    "Epilepsy": 'epilepsy[Title/Abstract]',
    "Fibromyalgia": 'fibromyalgia[Title/Abstract]',
}

# US patients (low, high), NIH funding in millions of dollars around 2017-2019: see docs/08 for the sources
PATIENTS = {"ME/CFS": (1_726_000, 3_746_000), "Multiple sclerosis": (486_000, 486_000)}
FUNDING_M = {"ME/CFS": 15, "Multiple sclerosis": 111}


def pubmed_count(term, year):
    r = requests.get(EUTILS, timeout=60, params={
        "db": "pubmed", "term": f"({term}) AND {year}[pdat]", "retmode": "json", "rettype": "count"})
    r.raise_for_status()
    return int(r.json()["esearchresult"]["count"])


def main():
    rows = []
    for label, term in TERMS.items():
        for y in YEARS:
            rows.append({"condition": label, "year": y, "papers": pubmed_count(term, y)})
            time.sleep(0.4)          # stay under the 3 requests per second limit
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "results" / "background_pubmed_counts.tsv", sep="\t", index=False)

    wide = df.pivot(index="year", columns="condition", values="papers")
    last5 = wide.loc[2021:2025].sum()
    per_patient = {}
    for c, (lo, hi) in PATIENTS.items():
        per_patient[c] = (last5[c] / (hi / 1000), last5[c] / (lo / 1000))
    funding = {c: (FUNDING_M[c] * 1e6 / PATIENTS[c][1], FUNDING_M[c] * 1e6 / PATIENTS[c][0]) for c in PATIENTS}

    summary = pd.DataFrame({
        "papers_2021_2025": last5[list(PATIENTS)],
        "papers_per_1000_patients_low": [per_patient[c][0] for c in PATIENTS],
        "papers_per_1000_patients_high": [per_patient[c][1] for c in PATIENTS],
        "nih_dollars_per_patient_low": [funding[c][0] for c in PATIENTS],
        "nih_dollars_per_patient_high": [funding[c][1] for c in PATIENTS],
    })
    summary.to_csv(OUT / "results" / "background_per_patient.tsv", sep="\t", float_format="%.2f")
    print(summary.round(2).to_string())

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for c in wide.columns:
        axes[0].plot(wide.index, wide[c], label=c, lw=2.2 if c == "ME/CFS" else 1.2)
    axes[0].set_ylabel("PubMed papers per year"); axes[0].legend(fontsize=7, frameon=False)
    axes[0].set_title("Research output")
    names = list(PATIENTS)
    lo = [funding[c][0] for c in names]
    hi = [funding[c][1] for c in names]
    axes[1].bar(names, hi, color=["#c0392b", "#7f8c8d"], alpha=0.35, label="fewer patients assumed")
    axes[1].bar(names, lo, color=["#c0392b", "#7f8c8d"], label="more patients assumed")
    axes[1].set_ylabel("NIH dollars per US patient per year")
    axes[1].set_title("Funding relative to the number of patients (about 2017-2019)")
    axes[1].legend(fontsize=7, frameon=False)
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUT / "figures" / "fig1_background.png", dpi=300, bbox_inches="tight")
    print(wide.loc[[2005, 2015, 2025]].to_string())


if __name__ == "__main__":
    main()
