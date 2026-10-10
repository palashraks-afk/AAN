# Is the possible discovery actually new? Deep check, 2026-10-10

Claim tested (DISCOVERY.md, A): "a calibrated, pre-registered test against 19 other traits shows that a small subset of brain cell types,
mainly amygdala excitatory neurons, carries ME/CFS signal beyond what generic brain traits show."
Method: about ten targeted searches this session plus page reads. Most sources are forum or blog posts read through a fetch/summary tool, so
every item below must be checked against the original before it is cited. Absence from search is not proof of absence.

## What already exists (so cannot be claimed as new)
| Item | Source | What it covers | Peer reviewed? |
|---|---|---|---|
| ME/CFS enriched in brain tissue (all 13 GTEx brain tissues of 54) | DecodeME preprint, medRxiv 2025.08.06.25333109 | tissue level | preprint |
| Neurons, not immune cells, in Siletti atlas; 24 of 3,172 region-by-type tests; central amygdala neurons in the top 10; eMSN claustrum | S4ME thread "Using MAGMA on ME/CFS genetic data" | cell types, Bonferroni, average-expression control | no |
| A 461-cluster MAGMA run (tralfamadorian97) described eMSN, amygdala excitatory and deep-layer intratelencephalic as separate signals | same thread (read through a summary; sources disagree on the exact wording; verify) | cell types | no |
| eMSN "best match", GTEx brain regions, comparison with other diseases from Duncan et al. 2025 (eMSN also appear for sleep duration, schizophrenia, depression, alcohol); post says the eMSN result is "far from unique" | ME/CFS Science blog, "Cell and tissue enrichment in ME/CFS" | descriptive comparison, no formal specificity test | no |
| Cell-type enrichment, pituitary tissue, glutamatergic synapses, rare-variant replication | Maccallini et al. 2026 (preprint) | tissue, cell type, gene sets | preprint |
| Heritability in inhibitory neurons in a fetal-derived atlas | Lee et al. 2026 (via S4ME) | different atlas | preprint |
| Same atlas-plus-MAGMA approach in other conditions (tinnitus: amygdala excitatory neurons; depression) | PMC12687797; PMC11829167 | other diseases | yes |

## What I did not find
A formal test of whether any ME/CFS cell-type signal exceeds what a panel of other traits shows (a specificity score against 19 traits), with
panel-robustness checks, conditioning on the panel's principal components, and a pre-registration. Duncan-style cross-disease plots are
descriptive. This is the only candidate for "new", and it is a method contribution, not a new cell type.

## How strong is the evidence for the amygdala result itself?
| Check | Result |
|---|---|
| Main test (model A) | Amex_153 z 3.4, FDR 0.033; Amex_175 z 3.2, FDR 0.040 |
| Stricter model B (also condition on neuronal expression) | no cluster passes FDR < 0.05 (best 0.067) |
| Panel-robust specificity (leave-one-out, half panels) | only the two amygdala clusters pass |
| My R1' (condition on 5 principal components of 19 other traits, MAGMA) | no cluster FDR < 0.05 (best 0.073); Amex_153 has the highest specificity score (3.69) |
| My GTEx check (v1.10, independent bulk data, different donors) | **does not replicate**: Brain - Amygdala is significant against all tissues (z 3.6, FDR 0.001) like every brain tissue, but among brain regions it is not special (p 0.49, specificity 0.25); pituitary p 7.6e-4 but specificity 1.1 (< 2); nerve p 0.75 |
| Second ME/CFS GWAS | same cases, so not independent |

## Verdict
- Not a discovery in the strong sense. The cell types are not new; some were named informally and by Maccallini-type analyses.
- What is defensible: a pre-registered, panel-calibrated, multi-method test, in which a small set of amygdala excitatory clusters keeps the
  highest specificity score, while the broader neuronal signal is shared with other brain traits. The result is modest (FDR 0.03 to 0.05), does not
  survive the stricter conditioning models, and is not seen in bulk GTEx amygdala (bulk tissue may hide a small cell population).
- Scores: novelty stays 3, becoming 4 only if the 17 Oct search is clean and the calibration is framed as the contribution. Do not write "discovery".
- Suggested wording: "to my knowledge, a panel-calibrated specificity test for ME/CFS cell types was not found in searches of 10 Oct 2026; the
  amygdala excitatory clusters keep the highest specificity score but the evidence is modest."
