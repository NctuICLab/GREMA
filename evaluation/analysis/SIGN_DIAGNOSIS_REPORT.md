# Why are the signs of GREMA wrong on real IRMA? A diagnosis with the true structure told

Earlier reports: the signs found on real IRMA are right in 25-50% of the called true edges (`SIGNS_REAL_DATA_REPORT.md`), not explained by the amount of data, the topology, a protein stage, an external input, a post-translational inhibition, or the Hill form (a linear model is as wrong; `NONHILL_FORM_REPORT.md`).
This report asks what is left: **when GREMA is told the true structure, does it find the right sign?** Scripts: `irma_oracle_knowledge.py`, `eval_irma_oracle.py`, `irma_sign_oracle_linear.py`. Archive: `evaluation/analysis/results/irma_sign_diagnosis_results.tgz` (see `evaluation/analysis/results/README.md`).
The gold standard is used on purpose: this is a diagnosis of why the sign is wrong, not a method. HFODE, `-f 1`, `-g 500`, 5 GA runs per setting.

**Answer: the wrong signs come from the data, not from the model.** With the true structure told, GREMA finds the right sign for only 43-50% of the signed edges, the sign of an edge is the same in every GA run, it is the sign of a plain linear regression on the same structure (14 of 14), and it changes between switch-on and switch-off for 4 of the 7 edges.
The mRNA time courses of the two IRMA experiments do not carry one consistent sign for these regulations.

## 1. D2: GREMA with the true structure (`irma_oracle_knowledge.py`, `eval_irma_oracle.py`)

The structure is the gold standard of IRMA (8 regulations); the signs of 7 are known (`input/IRMA/IRMA_signs.tsv`), G2 -> G4 has no known sign. Two kinds of knowledge files:
**one sign open**: the edge under test is `?`, the other six signed edges carry their true sign, G2 -> G4 is `?`, every other pair is 0; **no sign told**: the 8 true edges are `?`, every other pair is 0.
80 runs: 2 conditions x (7 edges + the no-sign setting) x 5 runs. Numbers: runs (of 5) with the right sign / wrong sign (no run left an edge uncalled).

| Edge (known sign) | switch-on, one sign open | switch-on, no sign told | switch-off, one sign open | switch-off, no sign told |
|---|---|---|---|---|
| G1 -> G2 (+) CBF1 -> GAL4 | 0 / 5 | 0 / 5 | 5 / 0 | 5 / 0 |
| G2 -> G3 (+) GAL4 -> SWI5 | 5 / 0 | 5 / 0 | 0 / 5 | 0 / 5 |
| G3 -> G5 (+) SWI5 -> ASH1 | 5 / 0 | 5 / 0 | 5 / 0 | 5 / 0 |
| G3 -> G1 (+) SWI5 -> CBF1 | **5 / 0** | 0 / 5 | 0 / 5 | 0 / 5 |
| G3 -> G4 (+) SWI5 -> GAL80 | 0 / 5 | 0 / 5 | 5 / 0 | 5 / 0 |
| G5 -> G1 (-) ASH1 -> CBF1 | 0 / 5 | 0 / 5 | 0 / 5 | 0 / 5 |
| G4 -> G2 (-) GAL80 -> GAL4 | 0 / 5 | 0 / 5 | 5 / 0 | 5 / 0 |

Over the 14 edge-condition cells (70 runs per setting): **one sign open 35 right, 35 wrong (50%)**; **no sign told 30 right, 40 wrong (43%)**; no run left an edge uncalled.

- **The sign is decided by the data, not by the search**: in all 28 cells (edge x condition x setting) all 5 GA runs give the same sign (the vote is 0.95-1.00); the genetic algorithm is not the source of the errors.
- **Telling the structure and the other signs does not fix the signs**: 50% is the level of chance. The one cell that changes with the knowledge is G3 -> G1 in switch-on (wrong when no sign is told, right when the others are told): there the sign of SWI5 -> CBF1 depends on whether ASH1 -> CBF1 is fixed.
- **The sign flips between the conditions** for 4 of the 7 edges (G1 -> G2, G2 -> G3, G3 -> G4, G4 -> G2), and G3 -> G5 is right in both; G3 -> G1 and G5 -> G1, the two edges into CBF1, are wrong in both conditions.
  So the slots (one regulator per gene when the structure is not told) and the competition between collinear regulators are not the explanation either: with the structure told they are gone.

## 2. D1: a linear regression with the true structure (`irma_sign_oracle_linear.py`, no GREMA)

OLS of the change `x_i(t+1) - x_i(t)` on a constant, `x_i(t)` and the true regulators (at lag 0, 1, 2), per condition. At lag 0 it gives the right sign for 2 of 7 edges in switch-on and 4 of 7 in switch-off (3 of 7 for both series with a constant for the condition; lags 1 and 2: 3-5 of 7).
ASH1 -> CBF1 (G5 -> G1) has the wrong sign at **every lag and in every condition**, with large t values (up to +3.2 in switch-on at lag 0 and +5.2 in switch-off at lag 1); G1 -> G2 in switch-on at lag 0 has t = -4.7.

**The regression and GREMA give the same sign for 14 of 14 edges and conditions** (lag 0, 5 of 5 GA runs each). The sign that GREMA reports on IRMA is the sign of the linear association of the one-step change with the regulator, given the structure; the Hill function adds nothing to it and removes nothing.

## 3. D3: the medium as an input (`INPUT_TERM_REPORT.md`)

Both conditions fitted together with the medium as an external input (5 runs, structure not told): signs right in 5 of 26 called true edges (19%) against 5 of 16 (31%) without the input. No improvement; G3 -> G5 right in all runs, G1 -> G2, G3 -> G1, G3 -> G4 and G5 -> G1 wrong in all.

## 4. What this says

- The signs of GREMA on this data set are those of a linear association; they are chance-level (43-50%) even with the true structure, and they depend on the condition. **The information about the sign of these regulations is not in a form that a one-step fit of the mRNA time courses of one experiment can read.** This is a property of the data (and of the sign table, see below), not of the Hill form, the slots, the optimiser or the number of runs.
- A flip of the sign between two conditions of the same network is an observable warning: **a sign that is not the same in the series of different conditions should not be reported as a regulation sign.** This could be a feature of GREMA (the sign stability across series, reported with the sign); it is not implemented.
- Candidate explanations that this analysis does not test: (a) the mRNA of a regulator is not the activity of the regulator (Swi5 and Ash1 are controlled at the protein level and by the cell cycle; the two edges into CBF1 are wrong in both conditions; GAL80 -> GAL4 acts on the protein); (b) the transition between the media dominates the co-movement of the mRNAs, so that associations reflect the common response and its timing and not the regulation;
  (c) the table of signs is from a secondary description of the network and may not describe what the mRNA does in these two experiments. Experiments that would separate them: perturbations of single regulators (the IRMA paper has them), the protein activity, or time courses with several different stimuli.

## Limits

- 7 signed edges on one network, 2 conditions of one series each, 5 GA runs per setting; at lag 0 the regression has 14-20 points and 3-4 coefficients, so single signs with small t values are weak (the large t values of ASH1 -> CBF1 are not).
- The one-step regression uses the mRNA at the same time or 1-2 steps earlier; the real sampling times (`-time`) were not used here.
- The gold standard is used in the knowledge files on purpose, so the figures are not a performance of the method.
