# One-page summary of my project

**Question.** Myalgic encephalomyelitis/chronic fatigue syndrome (ME/CFS) is a poorly understood illness with a recent large genetic study (DecodeME, 15,579 patients and 259,909 controls). Earlier work says the genetic risk sits in neurons. I asked a narrower question: which brain cell types carry ME/CFS risk that other brain traits do *not* show?

**Data.** All public: the DecodeME summary statistics, the Human Brain Cell Atlas (Siletti 2023, 461 cell clusters) and summary statistics for 19 comparison traits, including schizophrenia, Alzheimer's disease, depression and height. No patients or animals were involved in my own work.

**Methods.** I linked genetic signal to genes (MAGMA) and then to cell types, and built a specificity score that compares ME/CFS with the 19 other traits. I wrote the plan down before looking at results (tagged in my public repository), and every change is logged. Before reading the ME/CFS results, three control traits had to reproduce known biology (schizophrenia in neurons, Alzheimer's in microglia, height in connective tissue). I also planned to test the specificity score itself by using each comparison trait as the target in turn. I then ran robustness checks on the main result (other comparison panels, other definitions of specificity, a permutation test, and other settings) and checks in independent datasets.

**Results.** The three control traits reproduced the known biology, and the specificity score recovered the expected cell class in 4 of 5 known cases. For ME/CFS, 11 of 461 cell clusters were enriched, mostly excitatory neurons in the cortex and amygdala. Six of them were more specific to ME/CFS than the other traits, and the strongest was one amygdala cluster (Amex_153), which held up across every change to the comparison panel and the scoring definition. The same clusters also carry signal in chronic pain, though less strongly, and the signal was clearest in the infection-triggered patients.

**Limits of my study (not of earlier work).** The statistical evidence is modest. The signal weakens under stricter statistical models and with other gene-window settings, checks in independent datasets did not confirm the amygdala label in particular (each was a weak or indirect test), and the second DecodeME analysis shares the same patients, so it is not independent. It is genetic association, not causation, and a hypothesis for researchers, not a treatment finding.

**Relation to earlier work.** Neurons, and some of these cell types, were noted before in the DecodeME study, a 2026 preprint and informal analyses. What I could not find in the literature is a test of whether ME/CFS cell types are specific compared with a panel of other traits.

**Conclusion.** A small set of excitatory neurons, led by an amygdala cluster, carries ME/CFS genetic signal beyond what 19 other traits show. This gives researchers a more specific cell population to test than "the brain", and it offers a validated way to check whether other disease cell-type findings are specific.

Code, plans and results: https://github.com/palashraks-afk/AAN
