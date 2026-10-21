"""Build paper/BIBLIOGRAPHY.md from DOIs and PubMed ids (metadata from Europe PMC), check every link, and put the same list into paper/REPORT.md."""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
PAPER = Path("D:/AAN/paper")

DOIS = """10.1038/nature15393 10.1038/s41588-022-01024-z 10.1111/j.2517-6161.1995.tb02031.x 10.1038/s41588-020-0610-9
10.1038/s41586-018-0579-z 10.1038/ng.3211 10.1371/journal.pcbi.1004219 10.1101/2025.08.06.25333109 10.1038/s41588-018-0081-4
10.1038/s41593-018-0326-7 10.17226/19012 10.64898/2026.06.08.26355171 10.21203/rs.3.rs-9702020/v1 10.1073/pnas.1708274114
10.1038/s41588-021-00931-x 10.1126/science.add7046 10.1038/s41588-018-0129-5 10.1038/s41467-019-11181-1 10.1038/s41586-022-05275-y
10.2967/jnumed.113.131045 10.1038/s41586-021-03819-2 10.1093/nar/gkad1011 10.1038/s41467-025-63753-z 10.1016/j.neuron.2021.09.001
10.1126/science.aba7721 10.1126/science.aaz1776 10.1038/s41586-022-05473-8 10.1038/s41593-024-01834-w 10.1101/2021.11.24.469720
10.1016/s1474-4422(19)30320-5 10.1038/s41467-023-41210-z 10.1371/journal.pgen.1006125 10.1126/sciadv.abh2169
10.1093/nar/gkac1010 10.1093/nar/gkac1046 10.1093/nar/gkad1025""".split()
PMIDS = ["29892013", "33959723", "29942085", "34873335", "34741163", "35459240", "24390342", "28067908"]

MANUAL = [
    ("Benjamini Y, Hochberg Y.", 1995, "Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 57, 289-300", "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x"),
    ("Institute of Medicine.", 2015, "Beyond Myalgic Encephalomyelitis/Chronic Fatigue Syndrome: Redefining an Illness. Washington, DC: National Academies Press", "https://doi.org/10.17226/19012"),
    ("Brasher MS, et al.", 2023, "GWAS summary statistics for ICD-10 F41.9, anxiety disorder, unspecified. GWAS Catalog, study accession GCST90225527",
     "https://www.ebi.ac.uk/gwas/studies/GCST90225527"),
    ("Centers for Disease Control and Prevention.", 2024, "About ME/CFS. Accessed October 2026",
     "https://www.cdc.gov/me-cfs/about/index.html"),
    ("DecodeME Collaboration.", 2025, "DecodeME project and summary statistics, Open Science Framework. Accessed October 2026", "https://osf.io/rgqs3/"),
    ("Rakshit P.", 2026, "AAN: code, pre-registrations, results and figures for this project, GitHub repository", "https://github.com/palashraks-afk/AAN"),
]


def get(url, tries=3):
    for k in range(tries):
        try:
            return json.load(urllib.request.urlopen(url, timeout=40))
        except Exception as e:
            time.sleep(2 + k * 2)
    return None


def meta(query):
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=1&query=" + urllib.parse.quote(query)
    d = get(u)
    time.sleep(1.0)
    if not d or not d["resultList"]["result"]:
        return None
    return d["resultList"]["result"][0]


def fmt_authors(s):
    parts = [a.strip() for a in s.rstrip(".").split(",") if a.strip()]
    if len(parts) <= 4:
        return ", ".join(parts)
    return ", ".join(parts[:3]) + ", et al."


def fmt(r, doi):
    au = fmt_authors(r.get("authorString", ""))
    title = re.sub(r"<[^>]+>", "", r["title"]).rstrip(".")
    j = r.get("journalTitle") or ("Preprint" if doi.startswith(("10.1101", "10.64898", "10.21203")) else "")
    vol = r.get("journalVolume") or ""
    pg = r.get("pageInfo") or ""
    venue = f"*{j}*" if j else ""
    if vol:
        venue += f" {vol}"
    if pg and pg.lower() != "none":
        venue += f", {pg}"
    if doi.startswith("10.1101"):
        venue = "*bioRxiv/medRxiv preprint*"
    if doi.startswith("10.21203"):
        venue = "*Research Square preprint*"
    if doi.startswith("10.64898"):
        venue = "*medRxiv preprint*"
    if doi == "10.17226/19012":
        venue = "Washington, DC: National Academies Press"
    if doi == "10.1101/2025.08.06.25333109":
        au = "DecodeME Collaboration"
    if doi == "10.17226/19012":
        au = "Institute of Medicine"
    return au, r.get("pubYear"), title, venue.strip(), f"https://doi.org/{doi}"


entries, failed = [], []
for doi in DOIS:
    r = meta(f'DOI:"{doi}"')
    if not r:
        failed.append(doi)
        continue
    entries.append(fmt(r, doi))
for pm in PMIDS:
    r = meta(f"EXT_ID:{pm} AND SRC:MED")
    if not r or not r.get("doi"):
        failed.append("pmid " + pm)
        continue
    entries.append(fmt(r, r["doi"]))
print("could not resolve:", failed)

# link check through the DOI registry
bad = []
for e in entries:
    doi = e[4].replace("https://doi.org/", "")
    h = get("https://doi.org/api/handles/" + urllib.parse.quote(doi, safe="/"))
    if not h or h.get("responseCode") != 1:
        bad.append(doi)
    time.sleep(0.4)
print("DOIs not found in the registry:", bad)
for au, yr, ti, ve, url in [(m[0], m[1], m[2], "", m[3]) for m in MANUAL]:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        code = urllib.request.urlopen(req, timeout=40).status
    except Exception as ex:
        code = str(ex)[:40]
    print("web link", url, code)

lines = []
rows = [(e[0], e[1], e[2], e[3], e[4]) for e in entries] + [(m[0], m[1], m[2], "", m[3]) for m in MANUAL]
rows.sort(key=lambda x: re.sub(r"[^a-z]", "", x[0].lower()))
for i, (au, yr, ti, ve, url) in enumerate(rows, 1):
    body = f"{i}. {au} ({yr}). {ti}."
    if ve:
        body += f" {ve}."
    body += f" [{url}]({url})"
    lines.append(body)
text = "\n".join(lines)
(PAPER / "BIBLIOGRAPHY.md").write_text("# Bibliography\n\nReferences are listed alphabetically by first author. Every link was checked against the DOI registry or the web page on the date of submission.\n\n" + text + "\n", encoding="utf-8")

rep = (PAPER / "REPORT.md").read_text(encoding="utf-8")
i = rep.index("## References")
(PAPER / "REPORT.md").write_text(rep[:i] + "## References\n\n" + text + "\n", encoding="utf-8")
print(len(rows), "references written")
