# Pre-registration v1.16: MVP-only ME/CFS GWAS as a genuinely independent test (replaces the "descriptive meta-analysis" plan for M1 in v1.15 where possible)

Written 2026-10-10 (Mac session), tag `prereg-v1.16`. Written before the MVP-only file was downloaded or opened.

Why: the Zenodo meta-analysis file contains DecodeME and cannot test replication. The Million Veteran Program ME/CFS GWAS alone (PheCode 798.1, 3,891 cases, 439,202 controls, effective N 15,427; GWAS Catalog
GCST90479178) shares no participants with DecodeME or UK Biobank. It is underpowered and EHR-defined, and its sex mix differs (18.5% female cases), so a null is expected and weak evidence.

Rule (single test, same form as v1.11 H-R): in the MVP-only 461-cluster model A profile, the mean z of the six specific clusters (Amex_153, Amex_175, DLIT_152, DLIT_150, ULIT_121, Splat_402) exceeds the
95th percentile of 10,000 random six-cluster sets from the same profile (one-sided, p < 0.05). Also reported: how many of the six have z > 0 and the 11 main clusters' z.
Heritability gate: LDSC h2 z >= 4. Stated now: with 3,891 cases it will probably fail. If it fails the cluster profile is still computed but is labelled "below the heritability gate, descriptive only" and no
supported/not-supported call is made. If it passes the rule above decides.
The DecodeME + MVP meta-analysis (v1.15 M1) is still run as a descriptive comparison and is never called a replication.
Not allowed: changing the six clusters, the threshold or the permutation scheme after the MVP result is seen.
