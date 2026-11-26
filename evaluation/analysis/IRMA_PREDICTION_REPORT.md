# Can the model that GREMA fits on IRMA switch-on predict switch-off? (a fixed network, a new external input)

Question: GREMA gives a simulable ODE model. If it is right, it should predict the time course of the same network in a new condition. IRMA switch-on and switch-off are the same 5-gene network with a different medium (galactose / glucose).
Archive: `evaluation/analysis/results/irma_predict_results.tgz` (see `evaluation/analysis/results/README.md`). Script: `evaluation/analysis/scripts/predict_irma.py` (its rules are in the docstring and were fixed before any prediction was looked at).

**Answer: no. The fitted model fits its own condition well (skill 0.95 / 0.82) but does not predict the other condition, and it is not better than two trivial baselines in the reverse direction.**

## Method

- GREMA, HFODE, `-f 1 -g 500`, `IRMA_no_prior_knowledge.txt`, 5 GA runs per condition and per setting: the default number of slots, and `GREMA_SLOT_DIV=2` (two slots).
- Model of a gene: the run with the lowest fitness of its 30 runs (`out/G<i>.txt`: b, TransMax, n(1..5), k(1..5), degrade, fitness); `dx_i = TransMax_i prod_j 1/(1+(k_ij/x_j)^n_ij) + b_i - degrade_i x_i`, Euler step 1 per time point. The simulator reproduces the profiles that GREMA writes (`calprofile`) to 5e-7.
- **Free-running** prediction: all five genes simulated together from the observed first time point of the target series. **One-step** prediction: every gene from the observed values of the others, which is what GREMA fits; it is not a forecast and is given only for reference.
- Baselines: persistence (first value constant); training mean (the mean level of each gene in the training series, constant); no-regulation ODE (`dx = c - d x` per gene fitted on the training series, free-running).
- Metrics, mean over the five genes: RMSE; skill = 1 - MSE / MSE(persistence) (0 = as good as persistence, negative = worse); Pearson r of the predicted and observed course.

## Results (mean +- sd over the 5 GA runs)

| Train -> target | Method | RMSE | skill | r |
|---|---|---|---|---|
| on -> **off** | persistence | 0.286 | 0 | |
| | training mean | 0.160 | 0.29 | |
| | no-regulation ODE | **0.124** | **0.44** | 0.43 |
| | GREMA default, free-running | 0.605 +- 0.098 | **-4.69** +- 1.12 | -0.04 |
| | GREMA two slots, free-running | 1.326 +- 1.012 | **-53** +- 54 | -0.02 |
| | GREMA two slots, one-step | 0.152 | 0.14 | 0.45 |
| off -> **on** | persistence | 0.254 | 0 | |
| | training mean | 0.150 | 0.64 | |
| | no-regulation ODE | 0.136 | 0.64 | 0.42 |
| | GREMA default, free-running | 0.152 +- 0.002 | 0.41 | 0.33 |
| | GREMA two slots, free-running | 0.140 +- 0.005 | 0.55 +- 0.06 | 0.55 +- 0.04 |
| | GREMA two slots, one-step | 0.122 | 0.67 | 0.58 |
| on -> on (in-sample) | GREMA default / two slots, free-running | 0.048 / 0.042 | 0.95 / 0.96 | 0.89 / 0.91 |
| off -> off (in-sample) | GREMA default / two slots, free-running | 0.060 / 0.060 | 0.82 / 0.81 | 0.81 / 0.80 |

- on -> off: the free-running GREMA models are worse than keeping the first value (skill -4.7 and -53); in all 10 on -> off runs at least one gene goes above 1.5 (the observed maximum is 0.8; up to 11 in the two-slot runs), i.e. far outside the range of the data. Even the trivial no-regulation ODE (skill 0.44) is better.
- off -> on: the free-running skill (0.41 default, 0.55 two slots) is below the two trivial baselines (0.64); only the one-step prediction with two slots (0.67) is on a par with them, and it uses the observed values of the other genes.
- In-sample the same models fit well; the difference between the in-sample and the cross-condition skill is the part that the model does not extrapolate.

## What this does and does not show

- **The model has no input term.** Switch-on and switch-off differ in the medium, and the model only has gene-to-gene regulation; the two series start from different states and are driven differently. The test therefore shows that **the model as fitted by GREMA does not carry over to a new input**, not that a fixed regulation rule cannot. A model with an external input (a medium or a drug as `u(t)` in the production term) was not tried.
- The two series are one 5-gene system each, one series per condition, 5 GA runs; the baselines are weak by design. The no-regulation ODE and the training mean beat GREMA, which says the cross-condition difference is mostly in the levels, not the regulation.
- Not tested: forecasting inside one condition (fit the first part of a series, predict the rest). It separates "does not extrapolate in time" from "does not carry over to a new input"; GREMA runs on a prefix of the data are cheap.
- Consequence for drug perturbations: a model like this cannot be used to predict a new dose or drug unless the drug effect is expressed through measured genes (for instance as an initial value or an input to a measured node) and a held-out dose or drug is checked first as here.
