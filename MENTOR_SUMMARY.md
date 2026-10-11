# One-page summary of my project

**Question.** ME/CFS is a poorly understood illness with a recent large genetic study (DecodeME, 15,579 patients and 259,909 controls). Earlier work says the genetic risk sits in neurons. I asked a narrower question: which brain cell types carry ME/CFS risk that other brain traits do *not* show?

**Data.** All public: the DecodeME summary statistics, the Human Brain Cell Atlas (Siletti 2023, 461 cell clusters) and summary statistics for 19 comparison traits, including schizophrenia, Alzheimer's disease, depression and height. No patients or animals were involved in my own work.

**Method.** I linked genetic signal to genes (MAGMA) and then to cell types, and built a specificity score that compares ME/CFS with the 19 other traits. I wrote the plan down before looking at results (tagged in my public repository), and every change is logged. Three control traits had to reproduce known biology first: schizophrenia in neurons, Alzheimer's in microglia, height in connective tissue. They did. I also tested the specificity score by using each comparison trait as the target in turn; it recovered the expected cell class in 4 of 5 cases.

**Result.** 11 of 461 cell clusters are enriched for ME/CFS, mostly excitatory neurons in the cortex and amygdala. Six of them are more specific to ME/CFS than the other traits, and the strongest is one amygdala cluster (Amex_153). This holds up when I change the comparison panel, the scoring definition, and the number of genes it depends on. The signal appears mostly in the infection-triggered patients.

**Limits (I report these in the paper).** The evidence is modest. The signal weakens under stricter statistical models and with other gene-window settings, it was not replicated in two other datasets, and the second DecodeME analysis shares the same patients, so it is not independent. It is a hypothesis for researchers, not a treatment finding, and it is genetic association, not causation. The neuron and amygdala findings were noted informally before; what I could not find in the literature is this cross-trait specificity test.

Code, plans and results: https://github.com/palashraks-afk/AAN
