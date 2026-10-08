"""Write data/panel.tsv: the comparison traits used to calibrate the ME/CFS cell-type results.

Every trait is a harmonised GWAS Catalog summary-statistics file (GRCh38, rsIDs). The study table
and the harmonised file list come from EBI; both are downloaded by the first step below.
"""
import re
from pathlib import Path

import pandas as pd
import requests

PANEL = Path("D:/AAN_data/panel")
FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

# accession -> (short label, role in the analysis)
TRAITS = {
    "GCST90027158": ("alzheimer", "positive control (microglia)"),
    "GCST90245990": ("height", "negative control (non-neural)"),
    "GCST90018919": ("schizophrenia", "positive control (cortical neurons)"),
    "GCST005902": ("depression_broad", "brain trait / conditioning"),
    "GCST006476": ("neuroticism", "brain trait"),
    "GCST90018869": ("insomnia", "brain trait / conditioning"),
    "GCST90029007": ("bmi", "body size / conditioning"),
    "GCST90029013": ("education_years", "brain trait"),
    "GCST90038646": ("migraine", "brain trait"),
    "GCST009325": ("parkinson", "brain trait"),
    "GCST90027163": ("als", "brain trait"),
    "GCST90018840": ("epilepsy", "brain trait"),
    "GCST90275137": ("adhd", "brain trait"),
    "GCST90275151": ("autism", "brain trait"),
    "GCST003837": ("chronotype", "brain trait"),
    "GCST90225527": ("anxiety", "brain trait"),
    "GCST90077931": ("fibromyalgia", "related condition"),
    "GCST90080872": ("chronic_pain", "related condition"),
    "GCST90016564": ("ibs", "related condition"),
    "GCST90029070": ("crp", "immune/inflammation marker"),
    "GCST002318": ("rheumatoid_arthritis", "immune positive control"),
    "GCST90014023": ("type1_diabetes", "immune positive control"),
    "GCST004131": ("ibd", "immune positive control"),
}


def main():
    listing = requests.get(FTP + "harmonised_list.txt", timeout=180).text
    files = {}
    for line in listing.splitlines():
        m = re.search(r"/(GCST\d+)/harmonised/(.+\.h\.tsv\.gz)$", line)
        if m:
            files[m.group(1)] = line.strip().lstrip("./")

    studies = pd.read_pickle(PANEL / "harmonised_studies.pkl").set_index("STUDY ACCESSION")
    rows = []
    for acc, (label, role) in TRAITS.items():
        if acc not in files:
            print("no harmonised file for", acc)
            continue
        s = studies.loc[acc]
        s = s.iloc[0] if isinstance(s, pd.DataFrame) else s
        rows.append({
            "label": label, "accession": acc, "role": role, "n_total": int(s["N"]),
            "trait": s["DISEASE/TRAIT"], "first_author": s["FIRST AUTHOR"],
            "year": str(s["DATE"])[:4], "url": FTP + files[acc],
        })
    out = pd.DataFrame(rows)
    out.to_csv("D:/AAN/data/panel.tsv", sep="\t", index=False)
    print(out[["label", "accession", "n_total", "first_author", "year"]].to_string(index=False))


if __name__ == "__main__":
    main()
