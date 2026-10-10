# One-page summary for my mentor

**Project.** Which brain cell types carry the genetic risk for ME/CFS, and is that pattern specific to ME/CFS or shared with other brain conditions?

**Data.** DecodeME (UK, 15,579 patients, 259,909 controls, released 2025); the Human Brain Cell Atlas (Siletti 2023, 461 cell clusters);
summary statistics for 19 comparison traits. All public. No human or animal subjects were recruited by me.

**Method.** Gene-level association (MAGMA), cell-type enrichment, and a specificity score that compares ME/CFS with the 19 other traits. Planned
and written down before the results (pre-registration, with every change logged). Three control traits (schizophrenia, Alzheimer's, height) had to
reproduce known biology first, and they did.

**Result.** 11 cell clusters are enriched; six are more ME/CFS-specific than the other traits. Two amygdala excitatory clusters stay specific however
the comparison panel is changed. Several stricter checks weaken the result (no cluster passes when neuronal expression is also held constant), so it is
a hypothesis, not a confirmed discovery.

**What is not new.** Neurons, and some of these cell types, were reported informally before. The new part is the cross-trait specificity test.

**Limits.** Genetic association, not causation; three-donor atlas; European ancestry; the second analysis shares all patients.

All code, data notes, pre-registrations and results: https://github.com/palashraks-afk/AAN
