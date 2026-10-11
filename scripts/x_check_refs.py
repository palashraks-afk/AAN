"""M8: look up each reference by title/author/year in Crossref and Europe PMC and print the best match for manual checking."""
import requests

REFS = [
 ("Bulik-Sullivan 2015 Nat Genet", "LD Score regression distinguishes confounding from polygenicity in genome-wide association studies"),
 ("de Leeuw 2015 PLoS Comput Biol", "MAGMA: generalized gene-set analysis of GWAS data"),
 ("Siletti 2023 Science", "Transcriptomic diversity of cell types across the adult human brain"),
 ("Minikel 2024 Nature", "Refining the impact of genetic evidence on clinical success"),
 ("Watanabe 2019 Nat Commun", "Genetic mapping of cell type specificity for complex traits"),
 ("Bryois 2020 Nat Genet", "Genetic identification of cell types underlying brain complex traits yields insights into the etiology of Parkinson's disease"),
 ("Trubetskoy 2022 Nature", "Mapping genomic loci implicates genes and synaptic biology in schizophrenia"),
 ("Bellenguez 2022 Nat Genet", "New insights into the genetic etiology of Alzheimer's disease and related dementias"),
 ("Yengo 2022 Nature", "A saturated map of common genetic variants associated with human height"),
 ("Kurki 2023 Nature", "FinnGen provides genetic insights from a well-phenotyped isolated population"),
 ("Finucane 2018 Nat Genet", "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types"),
 ("Finucane 2015 Nat Genet", "Partitioning heritability by functional annotation using genome-wide association summary statistics"),
 ("DecodeME 2025 medRxiv", "Initial findings from the DecodeME genome-wide association study of myalgic encephalomyelitis/chronic fatigue syndrome"),
 ("Maccallini 2026 Research Square", "Biological Insights from Genome-Wide Association Studies and Whole Genome Sequencing of Myalgic Encephalomyelitis/Chronic Fatigue Syndrome"),
 ("Lee 2026 medRxiv", "Global and local genetic overlap among ME/CFS, irritable bowel syndrome and psychiatric traits: a hypothesis-generating analysis"),
 ("Slaughter 2026 medRxiv", "Seven replicated genomic associations of myalgic encephalomyelitis/chronic fatigue syndrome: a biobank study"),
 ("Cao 2020 Science", "A human cell atlas of fetal gene expression"),
 ("Tran 2021 Neuron", "Single-nucleus transcriptome analysis reveals cell-type-specific molecular signatures across reward circuitry in the human brain"),
 ("GTEx 2020 Science", "The GTEx Consortium atlas of genetic regulatory effects across human tissues"),
 ("Zhang 2022 Nat Genet scDRS", "Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data"),
 ("Fibromyalgia 2026 Nat Med", "The genetic architecture of fibromyalgia across 2.5 million individuals"),
 ("Taliun/TSPO PET 2014", "Evidence for widespread brain inflammation in chronic fatigue syndrome 11C-(R)-PK11195"),
 ("Institute of Medicine 2015", "Beyond Myalgic Encephalomyelitis/Chronic Fatigue Syndrome: Redefining an Illness"),
 ("Mirin Dimmock Jason 2020 Work", "Research update: The ME/CFS Illness"),
]
import re, time
for tag, title in REFS:
    q = re.sub(r"[^A-Za-z0-9 ]", " ", title)
    r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search", params={"query": f'TITLE:"{q[:140]}"', "format": "json", "pageSize": 5, "resultType": "lite"}, timeout=60).json()["resultList"]["result"]
    print("\n##", tag)
    print("   wanted:", title[:110])
    if not r:
        print("   NOT FOUND by exact title (check manually)")
    for it in sorted(r, key=lambda x: x["source"] != "MED")[:2]:
        print(f"   found : {it.get('title','')[:110]} | {it.get('authorString','')[:40]} | {it.get('pubYear')} | {it.get('journalTitle') or it['source']} | doi {it.get('doi')}")
    time.sleep(0.3)
