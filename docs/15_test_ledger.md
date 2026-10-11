# Test ledger: every pre-registered test and its outcome (multiple-testing accounting), 2026-10-10

Purpose: the headline FDR (0.03 to 0.05) is computed across the 461 clusters of one analysis; this table lists the separate pre-registered checks so a reader can see how many were run and how many came out supported.
"Supported" means the pre-registered rule was met. Tests marked R are robustness checks of the main result, I are independent or related-data checks, V are method validations. Details and numbers: `DISCOVERY.md`, `EXTENSION_LAYERS.md`.

| Layer | Pre-reg | Type | Tests | Outcome |
|---|---|:-:|---:|---|
| Control gate (schizophrenia, Alzheimer's, height) | v1 G3 | V | 3 | 3 of 3 passed (required before any ME/CFS result was read) |
| Localisation, model A | v1 T1 | main | 1 (461 clusters, FDR) | 11 clusters |
| Specificity against 19 traits | v1 T2 | main | 11 clusters | 6 specific |
| Conditioning | v1 T3 | R | 3 models | model B none at FDR < 0.05; model C 7 clusters |
| Confirmation (S-LDSC, permutation, locus/chromosome drop, second GWAS) | v1.3 | R | 4 per cluster | all 11 clusters confirmed (permutation p <= 0.0012 at 5,000 permutations); second GWAS shares cases |
| Method benchmark and validation of the score | v1.4, v1.12 | V | 4 methods + 5 classes | score recovers 4 of 5 expected classes |
| Groups, whole-body cell types, pathways, literature regions | v1.2, v1.3 | exploratory layers | 16 + 154 + 1,713 + 1 | none supported except two whole-body types (somatotrophs, one other) |
| Candidate-gene benchmark | v1 X3 | V | 3 diseases | 1 of 3 passed, so the gene list is exploratory only |
| Solidity checks W1, W3 to W9 | v1.14 | R | 8 | W1, W8, W9 pass; W3 (stable intervals) and W6 (region score) fail; W4 shows Amex_175 is mostly cortical; W5 and W7 descriptive |
| Specificity after pain added to the panel; Mac/Windows concordance | v1.18 | R | 2 | 5 of 6 clusters stay specific; Mac and Windows agree (Spearman 1.000) |
| Gene families | v1.7 F1 | I | 8 | 0 supported |
| Conditioning on principal components of 19 traits | v1.9 | I | 461 clusters | none at FDR < 0.05 |
| GTEx bulk tissue (amygdala, nerve, pituitary) | v1.10 | I | 3 | 0 supported |
| Independent amygdala/striatum atlas (2 regions, then 5) | v1.11, v1.17 | I | 2 + 2 | 0 supported |
| FinnGen malaise-and-fatigue cohort | v1.11 | I | 1 | not supported |
| Fetal whole-body atlas (peripheral groups, limbic neurons) | v1.13 | I | 6 + 1 | 0 supported |
| Replicated-loci genes | v1.15 | I | 2 | 0 supported |
| Related conditions (pain, fibromyalgia, IBS) | v1.15 | related | 3 | pain supported (shared pattern); the other two not |
| MVP-only ME/CFS GWAS | v1.16 | I | 1 | below the heritability gate; descriptive |

## What the ledger says
- Of the independent or related-data checks that test the ME/CFS hypothesis directly (about 22 individual pre-registered tests above), none was supported; the one supported test (pain) shows the pattern is shared, not replicated in ME/CFS.
- If one counts the main pipeline's 6 specific clusters as one finding, the number of separate tests run around it is large, so a single FDR of 0.03 to 0.05 should be read as modest. The report states this.
- The robustness checks mostly pass (W1, W8, W9, permutation, concordance), which supports that the result is a stable property of the data and the method, not an accident of one analysis choice. It does not make it replicated in independent data.
