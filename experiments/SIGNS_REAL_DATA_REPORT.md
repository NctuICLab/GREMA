# Are the signs of GREMA right on real data? Why not? A data check before the run

Three things: (1) the signs of the regulations found on the real IRMA data, for which the signs of 7 true regulations are known; (2) four checks of what could make them wrong; (3) `experiments/scripts/check_data.py`, a check of an input file for the properties that decide
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

- **On the real IRMA data the signs of GREMA are not better than random, and worse than "every edge `+`".** This contradicts the picture of the simulations (`SIGN_SIMULATION_REPORT.md`: 92-100% right when the model is right or mildly wrong).
  Section 2 tests what could cause it; none of the tests explains it.
- Consequence: the signs of GREMA are an advantage over the fast methods (which give no signs) in simulations only, not on real data.
- Limits: 7 edges, the signs come from a secondary description, 4 independent called edges per setting.

## 2. Why are the signs wrong? Four checks (archive `experiments/results/irma_sign_factors_results.tgz`)

**2a. The model is right, the topology and the amount of data are those of IRMA.** Data simulated with HFODE itself on the IRMA topology (`experiments/scripts/sim_gen.py <seed> <replicates> 0.05`: 8 edges with signs, 6 activating, 5 genes, 16 points, 5% noise, random initial values, 10 random parameter sets),
GREMA with the default settings (1 slot, `-f 1`, `-g 500`):

| Data | signed true edges called | sign right | vote AUROC (per run) |
|---|---|---|---|
| model right, **1 replicate** (like the real IRMA) | 37 / 80 | **35 / 37 (95%)** | 0.68 |
| model right, 5 replicates | 47 / 80 | 47 / 47 (100%) | 0.80 |
| real IRMA (`-f 1`, switch-on / switch-off) | 40 / 70 | 10 / 40 (25%) / 20 / 40 (50%) | 0.78 / 0.61 |

With the same topology and one series the signs are right (95%) when the data follow the model: the wrong signs on the real data are not caused by the amount of data or by the feedback loops of the topology.

**2b. An IRMA-like simulation with three switches** (`experiments/scripts/sim_irma_like.py`, `eval_irma_like.py`; the 5 genes and 7 signed regulations of IRMA, the mRNA is measured, 3 series of 20 points, 5% noise, 10 random parameter sets per condition, 8 conditions, GREMA default settings, 80 runs).
The switches: **protein stage** (the regulators act through a protein that follows the mRNA with a delay, rate 0.1-0.3), **external input** (a shift at step 5 that opens the regulation GAL4 -> SWI5: factor 0.2 + 0.8 u(t)), **post-translational inhibition** (GAL80 lowers the activity of the Gal4 protein instead of repressing the transcription of GAL4).
All the regulation functions are of the Hill type.

| Protein stage | Input | Post-transl. | AUROC | the 6 transcriptional edges called | sign right |
|---|---|---|---|---|---|
| off | off | off (the assumption of GREMA) | 0.819 | 44 / 60 | 44 (100%) |
| off | off | on | 0.809 | 43 / 60 | 43 (100%) |
| off | on | off | 0.802 | 39 / 60 | 39 (100%) |
| off | on | on | 0.808 | 42 / 60 | 42 (100%) |
| on | off | off | 0.665 | 29 / 60 | 29 (100%) |
| on | off | on | 0.556 | 20 / 60 | 20 (100%) |
| on | on | off | 0.680 | 28 / 60 | 28 (100%) |
| on | on | on | 0.590 | 20 / 60 | 20 (100%) |

The protein stage lowers the detection (AUROC 0.81 -> 0.62 on average, and fewer edges are called), the other two switches much less (input: 0.72 vs 0.71; post-translational inhibition: 0.69 vs 0.74). **None of them changes the sign** (265 of 265 called edges right).
GAL80 -> GAL4 (a transcriptional repression; a protein-level inhibition when the third switch is on) is the only regulation whose sign changes with a switch: with the third switch off it is found in 16 of 40 runs and always right, with it on in 6 of 40 runs and right in 2 of them.

**2c. A linear model on the real IRMA data** (`experiments/scripts/linear_signs_irma.py`; the change of a gene regressed on the other genes, ridge lambda 0.1 or one correlation at a time; the sign of the coefficient is the predicted sign): it gets the sign right for 2-3 of the 7 signed edges in switch-on and 3-4 in switch-off
(random 3.5, "every edge `+`" 5). **A model that has no Hill function is as wrong as GREMA on this data.** (The wrong edges are mostly the same: CBF1 -> GAL4, ASH1 -> CBF1, SWI5 -> CBF1.)

**2d. What this says about the Hill function.** The simulations of the sections 2a and 2b cannot show that the Hill function describes real regulation, because they were generated with Hill-type regulation functions; they show that GREMA recovers the signs when the regulation is of that type, also with a delay, an input and a protein-level inhibition.
The linear model on the real data shows that the wrong signs are **not specific to the Hill function**: a linear model is as wrong, which points to the data (confounding by a common regulator, e.g. SWI5 activates both ASH1 and CBF1, so ASH1 and CBF1 rise together although ASH1 represses CBF1; a few edges of one network; mRNA data for regulations that are partly at the protein level) rather than to the form of the model.
This confounding was not tested directly: in the simulations of 2a and 2b with the same topology the signs were right. **Not tested: regulation functions that are not Hill-type (non-monotone, a dual effect of the same regulator).**

## 3. A check of the data (`experiments/scripts/check_data.py`)

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
- The simulations of section 2 use Hill-type regulation functions (they cannot test whether the Hill function describes real regulation); 10 parameter sets per condition, 1 run each; the IRMA-like simulation has 3 series of 20 points (the real data: one series of 16 or 21).
- No sign check on the Hackett data (no signs in the gold standard) or on SOS (all true regulations are repressions: it is the case where the signs are known to be wrong).
