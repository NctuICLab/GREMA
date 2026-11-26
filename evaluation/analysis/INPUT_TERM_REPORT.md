# An external input for HFODE (`-input`): implementation and tests

Design: `evaluation/analysis/INPUT_TERM_DESIGN.md`. Motivation: the model that GREMA fits on one IRMA condition does not carry over to the other (`IRMA_PREDICTION_REPORT.md`, free-running skill -4.7 against persistence), and a stimulus that drives many genes is learnt as regulation between them (SOS).
Archive of inputs and results: `evaluation/analysis/results/input_term_results.tgz` (see `evaluation/analysis/results/README.md`). HFODE, `-f 1`, `-g 500`, default slots, 4 threads.

**Answer: the input term works as designed, and it finds the input and predicts a new condition when the data come from the model; on the two real data sets that were tried (IRMA, SOS) it gives no improvement.**
The criteria of the design were fixed before the runs; where a criterion failed it is reported as failed.

## 1. What was implemented

- `GREMA_main.pl -input U1,U2`: the named rows of the profile are external inputs. An input is a regulator whose time course is **given** (the same Hill factor, sign convention, slots, vote and `+ - ? 0` knowledge as a gene) and is never a target. The genes keep their numbers (the inputs are the last variables); `final_results.txt` has the lines `U1 G3 + ...`; the knowledge file has the inputs in the TF column only (an input without a line does not regulate).
  HFODE only, not with `-perturb`; `-time` is supported (the input is interpolated like the genes). **The solver is unchanged**: it already fits a gene from the observed values of all the other variables, so an input is a variable that is never a target.
- Without `-input` nothing changes: `check_input_regression.sh` compares the config, the knowledge files and the `run.log` mapping of the new and the previous version (identical). The runs themselves are random, so only the files are compared.
- Tools: `make_input.py` (a table of start time, level and baseline per series -> the input rows, scaled like the genes, `--shape decay`, `--reverse`), `check_data.py --input` (constant input, the same in every series, not on the scale of the genes, series of different length; a note instead of a warning when the series are similar because one stimulus drives them but the input differs between them),
  `grema_simulate.py` (the fitted model simulated for a new condition; it reproduces the profiles that GREMA writes to 5e-7 also with an input), the example `evaluation/analysis/examples/irma_input/`, and a README section.

## 2. Test 1: simulation with a known input (`sim_input_net.py`, `eval_input_sim.py`)

The IRMA topology (8 gene-gene regulations) plus a stimulus U1 that activates G1 and represses G4; 3 training series with the stimulus at (start, level) = (3, 1.0), (6, 0.6), (4, 0.8), one held-out condition (5, 0.7); 16 points, 5% noise, 10 parameter sets, one GREMA run each.

| | Result | Criterion |
|---|---|---|
| Input edges, AUROC (pooled, 10 seeds) | **0.948**; 15 of the 20 true input edges called, **15/15 signs right** | >= 0.8, sign >= 90%: **pass** |
| Gene-gene AUROC, with against without the input | 0.783 +- 0.079 against 0.763 +- 0.056 (+0.020) | not more than 0.05 lower: **pass** |
| Held-out condition, free-running skill (against persistence) | input **0.919 +- 0.123**, no input 0.601 +- 0.540, no-regulation ODE 0.428 +- 0.219, training mean -0.24; the input model beats the no-regulation ODE in 10 of 10 seeds (+0.49) | >= 0.1 above the no-regulation ODE: **pass** |
| Control: constant input | input-edge AUROC 0.453, 0 of 20 true edges called | near 0.5: **pass** |
| Control: input shuffled between the series | input-edge AUROC **0.702**, 6 of 20 called, 4 of 6 signs right | near 0.5 (0.4-0.6): **fail** |

The shuffled control fails because it is not an independent null: the three inputs are steps that rise at the points 3-6 and the shuffled input of a series correlates 0.73 with its true input (0.62-0.83), so part of the edge signal stays. The held-out prediction separates them clearly: the model trained on the shuffled input gets a skill of -5.3 (+- 18) against 0.92 with the true input.
The criterion is kept as failed; a control with an independent input (the stimulus permuted in time) was not run.
All the generators of this repository are of the Hill type, so this test shows that the implementation recovers a known input when the model is right, not that real regulation is of this form (compare the signs on real IRMA, `SIGNS_REAL_DATA_REPORT.md`).

## 3. Test 2: IRMA switch-on and switch-off as two series, the medium as input (`examples/irma_input/`)

The medium is constant within a series (1 in switch-on, 0.01 in switch-off): it is a label of the condition, so a model fitted on one condition sees one value of the input and **the held-out-condition test of `IRMA_PREDICTION_REPORT.md` cannot be repeated with an input** (a change of the design). Instead both series are fitted together and the later points are forecast (5 GA runs; the models free-running from the last training point, the input of the series given):

| Fitted on | Forecast | With the input (RMSE) | Without the input | No-regulation ODE | Training mean | Persistence |
|---|---|---|---|---|---|---|
| first 12 points | points 13-16 | 0.081 +- 0.003 | **0.077 +- 0.001** | 0.089 | 0.106 | 0.078 / 0.067 |
| first 8 points | points 9-16 (the transitions) | 0.113 +- 0.002 | 0.109 +- 0.001 | **0.101** | 0.113 | 0.116 / 0.123 |

- **Failed** (criterion: the input model below the model without the input and below the no-regulation ODE): in both forecasts the input model is slightly worse than the model without it, and with 8 points the plain no-regulation ODE is the best of all (skill +0.17 against -0.08 and -0.09 for the GREMA models). As in `IRMA_PREDICTION_REPORT.md`, a model with regulation does not beat a model without it on the real IRMA data.
- The 12-point forecast has little power: the last 4 points are a plateau (the persistence error 0.07-0.08 is only somewhat above the noise floor of 0.04-0.06 that two noisy points give with a noise sd of 0.03-0.05, estimated from the second differences), so no method can be told from keeping the last value.
- Descriptive: with the 16 points, U1 -> G3 is called `-` in all 5 runs (vote 1.00), U1 -> G1 in 1 of 5, the other genes almost never.

## 4. Test 3: SOS with the UV dose as input (`prep_sos_input.py`, `eval_sos_input.py`)

The exposure is at the first time point, so a step is constant in a series; the input is a **decay** of the dose (raw 5 and 20, scaled to 0.2575 and 1; half-time 60 min), to have an input that changes within a series. Main design: **leave one experiment out** (fit on three experiments, which have both doses; predict the fourth), 3 GA runs per fold and mode.

| Model | held-out skill, mean of the 4 folds | |
|---|---|---|
| no input | 0.782 | |
| decay of the dose | 0.787 | +0.005 above no input; above the no-regulation ODE by 0.011: both pass **by margins far below the differences between the folds (0.67-0.86)** |
| the same decay in the reverse order of time (control) | **0.811** | should be worse than the decay by >= 0.1: **fail** (it is better by 0.024) |
| training mean / no-regulation ODE | 0.784 / 0.776 | |

- All models predict the held-out experiment equally well (even the training mean, 0.78), because the four experiments respond alike (series similarity 0.62-0.68 in `check_data.py`): this test has little power, and the input model does not use the time shape of the input (the reversed one is as good).
- **Single-dose design** (fit on the two experiments of one dose, predict the other; `check_data.py` calls the input not identifiable: it is the same in both training series): fitted on UV 5 and predicting UV 20, the input model has a skill of **-0.03** against 0.73 without the input; the other direction 0.79 against 0.81. An input extrapolated beyond the levels seen in the training (0.26 to 1.0) makes the prediction worse.
- Edges and signs: LexA (G2) -> G3, G4, G5 and G8 are still mostly `+` (LexA is a repressor); the input U1 is called a regulator of almost every gene (vote 0.36-0.97, mostly `+`, always `-` for G7), so it absorbs the common time course without being specific.
- The sensitivity to the shape (half-times 30 and 120 min, the step) was **not run**: with the main result unchanged by the input, 24 more runs were not expected to change the conclusion.

## 5. What this says

- **The implementation is correct and useful where the data allow it**: the input and its sign are found (15/15), a new condition is predicted (0.92), a constant input is recognised as empty, `check_data.py` flags the inputs that cannot be identified, and nothing changes without `-input`.
- **On real data there is no evidence of benefit**: IRMA forecasts are not better with the input, and on SOS the held-out experiment is predicted equally well with or without it; the control with the reversed input is as good, and LexA signs do not improve. The two tests have limited power (a plateau, experiments that look alike), so this is "no evidence", not "evidence of no effect".
- **What would give a real test**: a designed experiment in which an input changes within a series and differs between series (two or more doses or start times, at least 8 points), with a held-out condition. The real data sets here do not have it: IRMA has a constant medium per series, SOS a single exposure at time 0.
- **Extrapolation of the input is unsafe**: the model with an input extrapolates badly to input levels outside the training range (SOS, -0.03 against 0.73); a prediction of a new dose or drug needs training data around it.
- Changes of the design that were made during the implementation (all before the runs they concern): `-perturb` is not supported; the IRMA held-out-condition test was replaced by forecasts of the later points; the SOS single-dose test was replaced by leave-one-experiment-out as the main design; `input_sensitivity.sh` is replaced by `make_input.py --shape`.

## Limits

- 10 simulated parameter sets with one GA run each; 3-5 GA runs on the real data; 2 series (IRMA) and 4 experiments (SOS); the generators are of the Hill type.
- The shape and the scale of the input are chosen by the user (the SOS decay with the half-time 60 min was fixed in advance and not varied).
- Not tested: S-system, `-time` with an input in a run, a delay between the input and the mRNA (phase 2 of the design), an input that is measured and noisy.
