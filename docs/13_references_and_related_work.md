# Reference check and related-work table (task M8, 2026-10-10)

Every reference named in `paper/DRAFT_report.md` and the project docs was looked up by title in Europe PMC (PubMed plus preprint servers) with `scripts/x_check_refs.py`; the DOI below is the
one the database returned. "Verified" means title, first author, year and journal agree with what the repository says. Items that failed the first lookup were fixed by author and year (notes in the last column).

## A. Methods and resources cited in the report
| Short name | Full reference | DOI | Verified | Note |
|---|---|---|:-:|---|
| Bulik-Sullivan 2015 | Bulik-Sullivan BK et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet 2015 | 10.1038/ng.3211 | yes | |
| Finucane 2015 | Finucane HK et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. Nat Genet 2015 | 10.1038/ng.3404 | yes | |
| Finucane 2018 | Finucane HK et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nat Genet 2018 | 10.1038/s41588-018-0081-4 | yes | LDSC-SEG |
| de Leeuw 2015 | de Leeuw CA et al. MAGMA: generalized gene-set analysis of GWAS data. PLoS Comput Biol 2015 | 10.1371/journal.pcbi.1004219 | yes | |
| Watanabe 2019 | Watanabe K et al. Genetic mapping of cell type specificity for complex traits. Nat Commun 2019 | 10.1038/s41467-019-11181-1 | yes | an author correction exists (2020) |
| Bryois 2020 | Bryois J et al. Genetic identification of cell types underlying brain complex traits yields insights into the etiology of Parkinson's disease. Nat Genet 2020 | 10.1038/s41588-020-0610-9 | yes | apostrophe broke the first lookup |
| Zhang 2022 (scDRS) | Zhang MJ et al. Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data. Nat Genet 2022 | 10.1038/s41588-022-01167-z | yes | |
| Siletti 2023 | Siletti K et al. Transcriptomic diversity of cell types across the adult human brain. Science 2023 | 10.1126/science.add7046 | yes | |
| Cao 2020 | Cao J et al. A human cell atlas of fetal gene expression. Science 2020 | 10.1126/science.aba7721 | yes | |
| Tran 2021 | Tran MN et al. Single-nucleus transcriptome analysis reveals cell-type-specific molecular signatures across reward circuitry in the human brain. Neuron 2021 | 10.1016/j.neuron.2021.09.001 | yes | |
| GTEx 2020 | GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 2020 | 10.1126/science.aaz1776 | yes | v8 median expression used |
| Trubetskoy 2022 | Trubetskoy V et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. Nature 2022 | 10.1038/s41586-022-04434-5 | yes | PGC3 |
| Bellenguez 2022 | Bellenguez C et al. New insights into the genetic etiology of Alzheimer's disease and related dementias. Nat Genet 2022 | 10.1038/s41588-022-01024-z | yes | confirmed by DOI |
| Yengo 2022 | Yengo L et al. A saturated map of common genetic variants associated with human height. Nature 2022 | 10.1038/s41586-022-05275-y | yes | |
| Kurki 2023 | Kurki MI et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature 2023 | 10.1038/s41586-022-05473-8 | yes | an author correction exists |
| Minikel 2024 | Minikel EV et al. Refining the impact of genetic evidence on clinical success. Nature 2024 | 10.1038/s41586-024-07316-0 | yes | the 2.6x figure |

## B. ME/CFS background cited
| Short name | Reference | DOI / source | Verified | Note |
|---|---|---|:-:|---|
| DecodeME 2025 | DecodeME Collaboration. Initial findings from the DecodeME genome-wide association study of ME/CFS. medRxiv 2025 | 10.1101/2025.08.06.25333109 | yes | preprint; no journal version found on 10 Oct |
| Maccallini 2026 | Maccallini P et al. Biological Insights from GWAS and WGS of ME/CFS. Research Square 2026 (full text and tables on Zenodo 20204356) | 10.21203/rs.3.rs-9702020/v1 | yes | not medRxiv; preprint |
| Lee 2026 | Lee JH. Global and local genetic overlap among ME/CFS, IBS and psychiatric traits. medRxiv 2026 | 10.64898/2026.06.08.26355171 | yes | preprint |
| Slaughter 2026 | Slaughter J et al. Seven replicated genomic associations of ME/CFS: a biobank study. medRxiv 2026 | 10.64898/2026.09.09.26362115 | yes | preprint |
| Kerrebijn 2026 | Kerrebijn I et al. The genetic architecture of fibromyalgia across 2.5 million individuals. Nat Med 2026 | 10.1038/s41591-026-04492-6 | yes | |
| Nakatomi 2014 | Nakatomi Y et al. Neuroinflammation in patients with CFS/ME: an 11C-(R)-PK11195 PET study. J Nucl Med 2014 | 10.2967/jnumed.113.131045 | yes | the draft wrote "JNM 55:945" |
| Mirin 2020 | Mirin AA, Dimmock ME, Jason LA. Research update: the relation between ME/CFS disease burden and research funding in the USA. Work 2020 | 10.3233/wor-203173 | yes | the repo title was paraphrased; use this one |
| IOM 2015 | Institute of Medicine. Beyond Myalgic Encephalomyelitis/Chronic Fatigue Syndrome: Redefining an Illness. National Academies Press, 2015 | NBK284897 | partly | the report itself is not in PubMed; cite the NAP page, not the book-review hit |

## C. Related work: what each covers and what it does not
| Paper | Data | Cell-type analysis? | Compared with other traits? | Conditioned on other traits? | Pre-registered? |
|---|---|---|:-:|:-:|:-:|
| DecodeME 2025 | DecodeME GWAS | no (13 brain tissues, 54 tested) | tissue contrast only | no | analysis plan on OSF |
| Maccallini 2026 | DecodeME + MVP (19,470 cases) | yes (FUMA, Siletti level 2 + white matter; eMSN, cerebellar glutamatergic) | no | average-expression control only | no |
| Lee 2026 | DecodeME, FinnGen IBS, PGC MDD, UKB loneliness | yes (Descartes fetal atlas: inhibitory, enteric neurons) | genetic correlations, not cell-type profiles | no | no |
| S4ME / ME/CFS Science analyses | DecodeME | yes (FUMA, 461-cluster MAGMA; eMSN, amygdala) | descriptive, other diseases from Duncan et al. 2025 | average expression | no |
| Kerrebijn 2026 | fibromyalgia, 54,629 cases | yes (neural cell types) | no | no | no |
| **This project** | DecodeME, 19 comparison GWAS, 4 atlases | yes (461 clusters) | yes (specificity score s) | yes (neuronal expression, PCs, three traits) | yes (v1 to v1.16) |
