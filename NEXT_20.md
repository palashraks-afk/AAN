# Next 20 things for the project side (set 2026-10-10)

Rules for all of them: write the test down (a short preregistration/PREREG_v1.x file, tagged) before running it; run it once; report whatever comes out, nulls included;
nothing here can add a new ME/CFS discovery claim, it can only make the existing one more solid, weaker, or better explained. Label anything not pre-registered as exploratory.
Keep result files in `results/`, write one short paragraph per task in `DISCOVERY.md` section 4 (or `EXTENSION_LAYERS.md` for the Mac), commit by file name, pull --rebase before pushing.

## Windows session (has all panel MAGMA output in D:\AAN_data; heavy compute and re-analysis)
W1. Power-matched calibration: scale or subsample the other traits' signal so every trait has the same cluster-level strength as ME/CFS, then recompute s. Answers "is ME/CFS specific only because it is weaker?"
W2. False-positive rate of the T2 rule: simulate null traits (same gene-level correlation structure, no cell-type signal) and count how often a cluster passes FDR < 0.05 and s >= 2.
W3. Chromosome-block bootstrap confidence intervals for s of the six specific clusters.
W4. One pre-specified hypothesis at dissection level: do amygdala excitatory clusters from the amygdala dissections carry the signal (nuclei named before looking), against the same clusters in other regions.
W5. Formal test of the six clusters in the infection-onset, non-infection, female and male subsets (24 tests, Bonferroni), to replace the descriptive counts.
W6. Amygdala region score for ME/CFS against the panel percentile (a single pre-set test, narrower than the failed literature region set).
W7. Power projection: with the final DecodeME case count (21,620 recruited vs 15,579 used), how large must the effect be for the amygdala clusters to pass model B and the PC conditioning? Gives the future-work statement.
W8. Alternative specificity definitions: rank-based s, top-decile specificity, trait-weighted panel (weights by profile similarity). Is the amygdala result the same?
W9. Driver analysis for Amex_153 and Amex_175: leave the top 1, 5 and 20 genes out, and leave each DecodeME locus out; does the signal rest on a few genes?
W10. Analysis-choice sensitivity: rerun the main test with different MAGMA SNP windows and gene annotations (0 kb, 35/10 kb, 50 kb) and report how s and the six change.

## Mac session (own downloads; independent data and the application package)
M1. Independent replication data: the Maccallini meta-analysis files on Zenodo (record 20204356, 312 MB) and any MVP-only ME/CFS summary statistics; test the six clusters in the part that does not contain DecodeME cases if it can be separated, or label it not independent.
M2. Independent amygdala atlas: the central amygdala single-nucleus atlas (Nature Communications 2026) or any open human amygdala snRNA-seq; does the Amex signature replicate?
M3. Replicated-loci genes (UK Biobank and All of Us biobank preprint: CLYBL, BICD1, GRIN2A, CSMD1, RORA and the others): are they more specific to the amygdala excitatory clusters than size-matched random gene sets?
M4. Related conditions: is the amygdala specificity also present in long COVID, fibromyalgia, chronic pain and IBS summary statistics (pre-registered, four tests)?
M5. Cross-atlas consistency: match the Siletti clusters Amex_153 and Amex_175 to another human brain atlas and check that the cell identity is the same.
M6. Reproducibility package: an environment file, a one-command rerun on a small test dataset, tests passing, a README for a stranger, and a data-availability statement ready for a DOI.
M7. Publication-quality figures: the calibration figure, the amygdala map, the validation figure and the controls in one clean multi-panel layout, with readable labels and colour-blind-safe colours.
M8. Reference check: verify every reference in `paper/DRAFT_report.md` against PubMed or the DOI, and build a related-work table (what each paper covers, what it does not).
M9. Judge Q&A pack: the 30 questions a neurologist on the panel is most likely to ask about this project, with short honest answers and the file that backs each one.
M10. Pre-registration timeline figure: planned versus changed versus deviations, from the git tags and `DEVIATIONS.md`, as one figure for the report (it shows the process was honest).
