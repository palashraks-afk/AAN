# Guardrails: patient safety, human-subjects question, AI use

## Patient-safety language (ME/CFS is a vulnerable community)
- Findings are hypotheses for researchers, not treatment recommendations. Never write "X could treat
  ME/CFS", dosing, or advice to use any drug. Say "genetic evidence supports prioritising X for
  laboratory/clinical investigation."
- Do not imply cure, diagnosis, or that genetics explains all ME/CFS. State h2 = 0.095 (modest).
- Cite the 2.6× genetic-support figure as an observational association (Minikel 2024), not a
  guarantee.
- Avoid language that psychologises or minimises the illness. Use "ME/CFS" and "myalgic
  encephalomyelitis" consistently; WHO ICD-10 G93.3 is neurological.
- If sharing with a patient organisation, send the report, not a "what you should take" list.

## Human-subjects question on the AAN form
The project analyses publicly released, de-identified, aggregate GWAS summary statistics and a public
cell atlas. No individual-level data, no recruitment, no contact with participants. Likely answer:
human subjects = No, animal subjects = No. **Confirm with the teacher/mentor and, if in doubt, ask
science@aan.com before answering.** Do not apply for individual-level DecodeME data (requires a
Data Access Committee proposal and is out of scope for the deadline).

## Data-use terms
Read the DecodeME OSF page terms; cite the DecodeME preprint and OSF project; do not attempt to
re-identify anyone. Cite the Human Brain Cell Atlas (Siletti et al. 2023, Science) and each GWAS used.

## AI-use disclosure (be honest; it protects the application)
- Claude assisted with idea search, literature scan, feasibility checking and the plan/repo scaffold.
- AAN requires original research and original written work. Palash runs the analyses, interprets them
  and writes the abstract, report and form answers himself.
- Permissible, disclosed: code debugging help, formatting. Not permissible: AI-written abstract/report
  prose submitted as own.
- Write one sentence in the report methods and one in the "who helped" answer.

## Pre-registration integrity
Freeze PREREG before downloading comparison-trait results. Commit and tag (`git tag prereg-v1`).
Any change after that is recorded as a dated deviation, never silently edited.

## Reproducibility
Public repo, pinned requirements, MD5-checked inputs, deterministic seeds, all excluded analyses
listed with reasons. Never claim a number that was not measured.
