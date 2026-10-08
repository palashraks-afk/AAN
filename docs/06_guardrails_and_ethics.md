# Guardrails and ethics

## Patient-safety language
ME/CFS patients have been burned by hype and by being dismissed, so the wording matters.

- Findings are hypotheses for researchers, not treatment advice. Never write "X could treat ME/CFS",
  never give doses, never tell anyone to take anything. Say "genetic evidence supports prioritising
  X for laboratory or clinical investigation."
- Don't imply a cure or a diagnosis, and don't imply genetics explains the whole illness. SNP
  heritability is modest (h2 = 0.095 on the liability scale).
- Quote the 2.6x genetic-support figure (Minikel et al., Nature 2024) as an association seen in past
  drug programs, not a promise for any single target.
- Use "ME/CFS" and "myalgic encephalomyelitis" consistently. WHO files it under neurological disease
  (ICD-10 G93.3).
- If the results go to a patient organisation, send the report, not a list of things to take.

## Human-subjects question on the AAN form
The project uses published, de-identified, aggregate GWAS summary statistics and a public cell
atlas. There is no individual-level data, no recruitment and no contact with participants, so the
likely answer is human subjects = No, animal subjects = No. Confirm with the teacher or mentor, and
if there is any doubt, ask science@aan.com before answering. Don't apply for DecodeME's
individual-level data; that needs a Data Access Committee proposal and won't fit the deadline.

## Data-use terms
Read the terms on the DecodeME OSF page and cite the preprint and the OSF project. Never try to
re-identify anyone. Cite Siletti et al. 2023 (Science) for the cell atlas and every GWAS used.

## Pre-registration
Freeze the thresholds in `preregistration/PREREG_v1.md` before looking at ME/CFS enrichment
results, then `git tag prereg-v1`. Any later change goes in `preregistration/DEVIATIONS.md` with a
date and a reason.

## Reproducibility
Public repo, pinned requirements, MD5-checked inputs, fixed random seeds, excluded analyses listed
with reasons. Don't report a number that wasn't measured.
