# Your to-do list (Palash)

Updated 2026-10-09. The computer does the analysis. These are the things only you can do, with dates. Deadline: **Tue 20 Oct,
11:59 p.m. CT**. Aim to submit **Mon 19 Oct**.

## Do these today (Fri 9 Oct): they depend on other people's time
- [ ] **Send the three signature requests** in the AAN portal (aan.secure-platform.com): parent or guardian, teacher (confirms
      you are a high-school student), mentor (confirms your work; if you have no mentor, a parent or teacher completes it). The
      system emails them. Without these the application cannot be processed.
- [ ] **Ask your teacher or mentor to confirm the "no human subjects, no animal subjects" answers.** The project uses published,
      de-identified summary statistics only. If in doubt, email science@aan.com before answering.
- [ ] **Keep the laptop plugged in with the lid open and sleep set to "never" while plugged in.** Jobs have been killed twice by
      sleep. Windows: Settings, System, Power.
- [ ] Decide who your **scientific mentor** is for the genetics part. Your oncology mentor may not be the right person for
      statistical genetics. Someone who has used GWAS or single-cell data can read your report and catch errors.

## Schedule
| Date | Computer (me) | You |
|---|---|---|
| Fri 9 Oct | height finishes, control check G3, first ME/CFS cell-type table | signature requests, mentor, power settings |
| Sat 10 Oct | other ME/CFS GWAS and panel running; confirmation checks | read `docs/08` background so you can explain the problem; look at the figures as they appear |
| Sun 11 Oct | discovery layers, evidence ledger, figures | read `DISCOVERY.md` and `docs/04`; list questions you cannot answer yet |
| Mon 12 Oct | panel finishes; decision gate: is the translational layer worth keeping? | start writing the Introduction in your own words (`paper/DRAFT_report.md` section 1 has the points to cover) |
| Tue 13 Oct | final numbers frozen in `results/` | write Methods in your own words (the draft describes what the code does; rewrite it) |
| Wed 14 Oct | | write Results from `results/EVIDENCE_LEDGER.md` and the figures; mentor reads a first draft |
| Thu 15 Oct | | write Interpretation and Limitations; fix mentor comments |
| Fri 16 Oct | | abstract (under 300 words), bibliography (check every citation), figure legends |
| Sat 17 Oct | **prior-art re-search** (terms in `HANDOFF.md` section 7); update `docs/02` | revise claims that changed; teacher reads the whole thing |
| Sun 18 Oct | | export abstract, report, bibliography as PDF; read them out loud once; check word count |
| Mon 19 Oct | | **submit.** Do the project-information answers (role, interest, who helped, funding, future, family, other interests) |
| Tue 20 Oct | backup only | deadline 11:59 p.m. CT; do not plan to use this day |

## Useful additions from the extension work (2026-10-10)
- In the report's Limitations, add one honest sentence from `docs/12_discovery_novelty_check.md`: the amygdala signal is modest, does not survive the
  stricter conditioning models, and does not replicate in bulk GTEx (which may hide a small cell population).
- Read the originals before citing: Maccallini et al. 2026, Duncan et al. 2025, and the S4ME / ME/CFS Science posts (all were read through summaries).
- 17 Oct prior-art search: add "amygdala", "calibrated", "specificity" and "panel of traits" to the terms; update `docs/12`.

## What you write yourself (the prize needs your own words)
- The abstract, the report and the project-information answers in the form. Drafts in `paper/` are scaffolding only.
- In "who helped", name every mentor and helper accurately.
- Be ready to explain every sentence: what MAGMA does, what a cell-type enrichment means, why the controls matter, what the
  pre-registration is, and what the result does *not* show.

## Words to avoid in the application
"first", "never done before", "groundbreaking", "cure", "treatment for ME/CFS", "proves". Use "to my knowledge not previously
reported (searched [date])", "suggests", "genetic evidence supports prioritising X for study".

## Questions a judge could ask (prepare answers)
1. Why did you pre-register, and what did you do when the DecodeME files turned out not to be quality-controlled? (`DEVIATIONS.md`)
2. How do you know the pipeline works? (controls: schizophrenia, Alzheimer's, height; `docs/10_self_review.md`)
3. Neurons are enriched for almost every brain trait, so why is yours interesting? (the specificity test against 23 traits)
4. What would make your result wrong? (association not causation; three-donor atlas; central nervous system only; European
   ancestry; samples overlap between the DecodeME analyses)
5. Could this help patients? (it can point researchers at cell types and genes; it is not a treatment and says nothing about
   what any patient should take)

## If something goes wrong
| Problem | Do this |
|---|---|
| Jobs stopped (sleep, restart) | `cd D:\AAN\scripts` then `python run_everything.py`; it skips finished work and carries on |
| Control check says FAIL | stop; debug the pipeline; do not read the ME/CFS result; note it in `docs/10_self_review.md` |
| Behind schedule on 13 Oct | drop the translational layer and the neglected-conditions extension; the core is T1, T2, the controls and the confirmation checks |
| Someone publishes the same analysis | cite it, say yours is an independent replication, lower the novelty claim in `SCORES.md` |
| A git commit looks wrong | `git status` and `git diff --stat` first; never `git add -A` blindly |
