"""Pull the targets of approved drugs for the benchmark diseases from the Open Targets Platform API.

Written to data/benchmark_targets.tsv (small, committed). Run it before ranking anything, so the
list of "known winners" can't be shaped by the results.

    python scripts/get_known_targets.py
"""
import pandas as pd
import requests

API = "https://api.platform.opentargets.org/api/v4/graphql"
# IBD is queried through its parent and its two main children: the parent node alone lists only a
# handful of approved drugs because most approvals are recorded under Crohn's disease or ulcerative colitis
DISEASES = {
    "migraine": ["MONDO_0005277"],
    "rheumatoid_arthritis": ["MONDO_0008383"],
    "ibd": ["MONDO_0005265", "MONDO_0005011", "MONDO_0005101"],
}

QUERY = """
query($efo: String!) {
  disease(efoId: $efo) {
    name
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          name
          mechanismsOfAction { rows { mechanismOfAction targets { id approvedSymbol } } }
        }
      }
    }
  }
}
"""


def main():
    out = []
    for label, efos in DISEASES.items():
        for efo in efos:
            r = requests.post(API, json={"query": QUERY, "variables": {"efo": efo}}, timeout=180)
            r.raise_for_status()
            disease = r.json()["data"]["disease"]
            rows = disease["drugAndClinicalCandidates"]["rows"]
            print(label, efo, disease["name"], "| candidates:", len(rows))
            for row in rows:
                if str(row["maxClinicalStage"]).upper() != "APPROVAL":
                    continue
                moas = (row["drug"].get("mechanismsOfAction") or {}).get("rows") or []
                targets = {t["approvedSymbol"]: (t["id"], moa["mechanismOfAction"])
                           for moa in moas for t in (moa["targets"] or [])}
                for symbol, (ensembl, mechanism) in targets.items():
                    out.append({"disease": label, "efo": efo, "target": symbol, "ensembl": ensembl,
                                "drug": row["drug"]["name"], "mechanism": mechanism,
                                "n_targets_of_drug": len(targets)})
    df = pd.DataFrame(out).drop_duplicates(["disease", "target", "drug"])
    df.to_csv("D:/AAN/data/benchmark_targets.tsv", sep="\t", index=False)
    selective = df[df["n_targets_of_drug"] <= 3]
    print("all approved-drug targets:\n", df.groupby("disease")["target"].nunique())
    print("targets of selective drugs (<= 3 targets):\n", selective.groupby("disease")["target"].nunique())


if __name__ == "__main__":
    main()
