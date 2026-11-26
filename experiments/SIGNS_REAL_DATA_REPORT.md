# Are the signs of GREMA right on real data? A data check before the run

Two things: (1) the signs of the regulations found on the real IRMA data, for which the signs of 7 true regulations are known; (2) `experiments/scripts/check_data.py`, a check of an input file for the properties that decide
whether a regulation can be inferred (used in the README section "Is my data suitable?").

## 1. Signs on IRMA (`experiments/scripts/eval_signs_real.py`, `input/IRMA/IRMA_signs.tsv`)

Signs of the true regulations (5 genes: G1 CBF1, G2 GAL4, G3 SWI5, G4 GAL80, G5 ASH1), from the model of the IRMA network of Cantone et al. (2009) in
the arXiv paper [1005.2107](https://arxiv.org/pdf/1005.2107) ("Efficient parameter search for qualitative models of regulatory networks using symbolic model checking"; its figure 1 and equations; the source is a description of the network, not the figure of the original paper):

| Regulation | sign | | Regulation | sign |
|---|---|---|---|---|
| CBF1 -> GAL4 | + | | SWI5 -> GAL80 | + |
| GAL4 -> SWI5 | + | | ASH1 -> CBF1 | - |
| SWI5 -> ASH1 | + | | GAL80 -> GAL4 (inhibition of Gal4 by Gal80) | - |
| SWI5 -> CBF1 | + | | | |

The 8th regulation of the gold standard, GAL4 -> GAL80, is not in that model; its sign is not used. So 7 signed edges, 5 activating: "every edge `+`" is right for 71% of them, random 50%.

GREMA (HFODE, `-f 1` / `-f 3`, `-g 500`, 10 runs per setting; `final_results.tgz` and `ssystem_slots_results.tgz` of `experiments/results/`):

| Setting | signed true edges called a regulation | of those, sign right |
|---|---|---|
| `-f 1`, switch-on | 40 / 70 | **10 / 40 (25%)** |
| `-f 1`, switch-off | 40 / 70 | 20 / 40 (50%) |
| `-f 3`, switch-on | 49 / 70 | 19 / 49 (39%) |
| `-f 3`, switch-off | 40 / 70 | 20 / 40 (50%) |
| `-f 1`, 2 slots, switch-on | 59 / 70 | 17 / 59 (29%) |
| `-f 1`, 2 slots, switch-off | 37 / 70 | 16 / 37 (43%) |

The sign of an edge is almost the same in all 10 runs of a setting (an edge is 10/10 right or 0/10 right), so each setting has about 4 independent edges, not 40. Systematic errors: **ASH1 -> CBF1** (true `-`) is called in every setting and always `+`;
**CBF1 -> GAL4** (true `+`) is called in 4 of the 6 settings and always with the wrong sign `-`; **SWI5 -> ASH1** (true `+`) is right in all settings (except 1 run of 10 in one).

- **On the real IRMA data the signs of GREMA are not better than random, and worse than "every edge `+`".** This contradicts the picture of the simulations (`SIGN_SIMULATION_REPORT.md`: 92-100% right when the model is right or mildly wrong): the model of the real data is far from the HFODE form.
  Possible reasons (not tested): IRMA measures mRNA, while the Gal4-Gal80 regulation is a protein-protein interaction and Cbf1 acts through a promoter that depends on the medium; one series of 16 / 21 points.
- Consequence: the signs and the dynamics of GREMA are not an advantage over the fast methods (which give no signs) that was shown on real data; this was shown in simulations only.
- Limits: 7 edges, the signs come from a secondary description, 4 independent called edges per setting.

## 2. A check of the data (`experiments/scripts/check_data.py`)

Prints the size of the data (series, points, samples = series x (points - 1)), the mean |correlation| between the genes and the **series similarity** (mean correlation of the same gene between two different series), and warns (rules of thumb, not tests):
fewer than 3 series, fewer than 8 points, fewer than 30 samples, mean |correlation| >= 0.55, series similarity >= 0.5.

| Data | series similarity | mean \\|correlation\\| | warnings |
|---|---|---|---|
| IRMA (switch-on) | one series | 0.37 | 1 series, 15 samples |
| SOS, 4 experiments | **0.66** | 0.49 | same time course |
| Hackett (first TF subset, 10 experiments) | 0.03 | 0.30 | 7 time points |
| DREAM4 size 10 #1 | 0.16 | 0.31 | none |
| simulation, random initial values | 0.10-0.23 | 0.39 | |
| simulation, shared driver (transcription / cleavage) | 0.96 / 0.80 | 0.56 / 0.64 | collinear, same time course |

SOS, the data on which both methods score well without using the regulation (`METHODS_REPORT.md` section 1), is the one real data set that is flagged; the simulations with a shared driver are flagged; data with different perturbations are not.
The thresholds were set on these data sets (SOS against Hackett, DREAM4 and the simulations) and not validated on others. A first version used the R^2 of a series predicted by the mean of the others; it did not flag SOS (0.38) and gave meaningless values (-155 on Hackett) and was replaced by the correlation.
The README section also repeats the other findings of the reports: the permutation control (`permute_replicates.py`), the number of series, the measured variable that should be the regulator's activity.

## Limits

- The sign check has 7 edges on one network (two switch conditions) and a secondary source for the signs; the check of the data was developed on the same data sets that motivated it.
- No sign check on the Hackett data (no signs in the gold standard) or on SOS (all true regulations are repressions: it is the case where the signs are known to be wrong).
