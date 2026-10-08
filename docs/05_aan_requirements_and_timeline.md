# AAN Neuroscience Research Prize 2027 — requirements and schedule

Sponsor: American Academy of Neurology + Child Neurology Society.
Application portal: https://aan.secure-platform.com/a/solicitations/682/home
Contact: Katie Anderson, Scientific Education Program Manager — science@aan.com

## Dates
- Applications opened 29 Jul 2026. **Close 20 Oct 2026, 11:59 p.m. CT.**
- Round-1 advancers notified late Nov 2026; finalists by end of Jan 2027.
- Prizes: three AAN winners $1,000 each + poster at the 2027 AAN Annual Meeting (Washington, DC),
  travel for winner + parent/mentor; one Bhuwan Garg High School Prize winner ($1,000, CNS meeting).
  Travel depends on the in-person meeting happening.

## Eligibility
Grades 9–12 in a US secondary school; original research AND original written work of the applicant;
individual project (no group projects); projects need not be in a formal lab; family of judges /
AAN Science Committee / AAN staff ineligible.

## Required materials
1. Completed application form (profile + project information).
2. **Abstract, ≤ 300 words** (.doc or .pdf). Round 1 is judged on this alone.
3. **Research report** (.doc or .pdf).
4. **Bibliography** (.doc or .pdf).
5. E-signature confirmation statements, requested through the portal:
   parent/guardian, teacher (confirms high-school student), mentor (confirms work; if no mentor,
   parent/guardian or teacher completes it). The system emails them — **send these requests first**.
6. Human / animal subject questions (yes/no each; forms required if yes).

## Judging
Round 1 (abstract only): relevance to neuroscience (behavior/psychology discouraged unless clear
neurophysiology link); abstract writing quality (clarity, methods, succinct); research potential
(innovation, impact, fills a gap).
Final round (full project): relevance; creativity (originality of problem-solving approach);
interpretation of data (feasible hypothesis/method, significance in perspective, **pitfalls
addressed**); research report (organised, well written, figures/tables labelled and readable).

### Mapping to this project
| Criterion | Where it is met |
|---|---|
| Relevance to neuroscience | Which cells of the human nervous system carry ME/CFS genetic signal; WHO-neurological disease |
| Abstract quality | 300-word structure: gap → question → method (pre-registered, calibrated) → result → limits |
| Research potential | Neglected, large population; hypothesis list for researchers |
| Creativity | Specificity calibration + benchmark-validated translational dossier |
| Interpretation / pitfalls | PLAN §7 limitations; sample overlap; null result pre-committed |
| Report quality | Figure plan below |

## Form fields (answers must be Palash's own words — prompts only)
Project title; role and how long on project; how you became interested and what you gained; who helped
and how (name mentors and anyone else who helped, accurately); financial support (none, if true); future plans;
family in science/medicine (yes/no); other areas of science that interest you; human subjects (see
docs/06); animal subjects (No). Profile fields are pre-filled from the AAN member profile.
Do not put phone, address or other personal data in this public repository.

## Day-by-day schedule (today = Wed 7 Oct)
| Date | Task |
|---|---|
| Wed 7 Oct | Repo + plan done. **Send the 3 signature requests** (need names/emails of parent, teacher, mentor). Choose mentor for genetics |
| Thu 8 Oct | Set up D:\AAN_data; install tools (MAGMA, LDSC-py3, R or Python coloc); re-run prior-art search; **freeze PREREG v1 and git-tag it**; download DecodeME files (after approval) and verify MD5 |
| Fri 9 Oct | Read README_shared_sum_stats.txt; QC + harmonise; LDSC h2 for every GWAS; resolve build/LD route |
| Sat 10 Oct | Build specificity matrix from agg.loom; assemble comparison-trait panel (download, QC); pipeline on positive controls (SCZ, AD, height) |
| Sun 11 Oct | MAGMA gene-level for ME/CFS and panel |
| Mon 12 Oct | Cell-type enrichment ME/CFS + panel; specificity calibration; LDSC-SEG confirmation |
| Tue 13 Oct | Conditioning (depression/BMI/insomnia); subtype and sex contrasts. **Go/no-go gate**: if core enrichment is not working, drop translational + stretch |
| Wed 14 Oct | Translational: coloc / MR inputs, druggability tables, safety flags |
| Thu 15 Oct | Benchmark on positive-control diseases; freeze all analyses |
| Fri 16 Oct | Figures + tables; draft report (Palash writes) |
| Sat 17 Oct | Abstract ≤ 300 words; bibliography; second prior-art re-search; mentor/teacher read |
| Sun 18 Oct | Revise; export PDFs; check word counts; test upload |
| Mon 19 Oct | Buffer; **submit** |
| Tue 20 Oct | Hard deadline 11:59 p.m. CT — do not plan to use this day |

## Report structure and figure plan
Abstract; Introduction (ME/CFS burden, neglect, gap); Methods (data, QC, gates, enrichment,
calibration, conditioning, dossier, benchmark, pre-registration); Results; Interpretation;
Pitfalls and limitations; Future work and how findings could be tested; References.
Figures: F1 design flowchart; F2 raw enrichment by supercluster (ME/CFS vs controls); F3 calibration
(ME/CFS percentile vs panel traits); F4 conditioning effect; F5 subtype/sex contrasts; F6 benchmark
(known drug targets rank); F7 evidence-tier table of candidate targets. Tables: data inventory,
panel traits, pre-registered thresholds, excluded analyses with reasons.

## Submission checklist
- [ ] Parent/guardian, teacher, mentor requests sent and completed
- [ ] Abstract ≤ 300 words (count it), PDF
- [ ] Report PDF with readable figures
- [ ] Bibliography PDF
- [ ] Human subjects / animal subjects answers confirmed with teacher
- [ ] "Who helped" answer lists every mentor and helper accurately
- [ ] No claims of "first/never/cure"
- [ ] Repo public, tagged, reproducible
