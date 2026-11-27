# Can augmented data take the place of real series? A simulation with a known truth

Question: when the data are not enough (one series), can the data augmentation of GREMA (`-perturb`, Gaussian-process samples) make up for it? On real IRMA it did not help (`IRMA_REPORT.md`, `input/IRMA/README.md`), but there the truth is not known.
Here the truth is known, and the augmented one-series data are compared with real series of the same systems.
Archive: `evaluation/analysis/results/sim_augment_results.tgz` (see `evaluation/analysis/results/README.md`). Scripts: `sim_aug_prep.py`, `eval_augment.py`, `gp_v2.py`. HFODE, `-f 1`, `-g 500`, default slots, no prior knowledge.

**Answer: no. None of the augmentations came close to a second real series, and the largest gain of an augmentation (+0.04 AUROC, not clear) is a quarter of the gain of one more real series (+0.10).**
The rules were fixed before the runs; all three augmentations fail them.

## Design

- 5 random signed networks of 10 genes (14-18 regulations each, `sim_signed_net.py`: the generating model is the model of GREMA, 3% noise, random initial states), 5 series of 20 points each; series 1-4 are the "real" data, series 5 is never used for a fit and tests the prediction.
- Conditions (per network): **1, 2, 3, 4 real series**; **1 series + `-perturb 3`** (4 series, 3 noisy copies; the noise 0.03 is the noise of the data) and **`-perturb 9`** (10 series); **1 series + 10 GP samples** (`gp_v2.py`, the better of the two GP versions of the IRMA analysis, hyperparameters marginalised, with observation noise; GREMA is run on every sample alone and an edge is scored by the fraction of the 10 runs that call it).
  80 GREMA runs. The augmentations are made from series 1 only; no gold standard and no other series is used to make them.
- Measures: AUROC of the vote score against the true edges of the 90 non-self pairs; the free-running prediction of the held-out series 5 (skill = 1 - MSE / MSE of persistence) and whether it is stable (skill > 0).
- Rule: an augmentation takes the place of real series if it closes at least half of the AUROC gap between 1 and 4 real series (paired over the networks) **and** its prediction is stable in at least as many networks as with 4 real series minus one.

## Results (5 networks)

| Condition | AUROC (mean +- sd) | paired gain over 1 real series (+- s.e.) | held-out skill, median | stable prediction |
|---|---|---|---|---|
| 1 real series | 0.759 +- 0.074 | | 0.63 | 4 of 5 |
| 1 series + `-perturb 3` (4 series) | 0.780 +- 0.075 | +0.020 +- 0.006 (5 of 5 networks) | 0.23 | 3 of 5 |
| 1 series + `-perturb 9` (10 series) | 0.757 +- 0.066 | -0.002 +- 0.023 | 0.70 | 4 of 5 |
| 1 series + 10 GP samples | 0.797 +- 0.103 (mean vote 0.788) | +0.037 +- 0.028 (3 of 5) | -0.37 | 2 of 5 |
| 2 real series | 0.856 +- 0.069 | **+0.096 +- 0.016** (5 of 5) | 0.83 | 4 of 5 |
| 3 real series | 0.889 +- 0.073 | +0.129 +- 0.021 | 0.92 | 5 of 5 |
| 4 real series | 0.901 +- 0.062 | **+0.141 +- 0.014** (5 of 5) | 0.90 | 5 of 5 |

By the rule: `-perturb 3` closes 0.020 of the 0.141 gap and has a stable prediction in 3 of 5 networks (4 real: 5): **does not take the place**; `-perturb 9` closes none (-0.002; prediction stable in 4 of 5): **does not**; the GP samples close 0.037 and are stable in 2 of 5: **does not**.

- **One more real series is worth +0.10 AUROC**, and the third and fourth add +0.03 and +0.01 more: the first additional series (a new initial state of the system) carries most of the information.
- **The built-in augmentation**: `-perturb 3` gives a small, consistent gain (+0.020 +- 0.006, in all 5 networks) but its held-out prediction is worse than with the single series (median skill 0.23 against 0.63); `-perturb 9` gives no gain although it has 10 series. A gain that does not grow with more copies is not new information (the copies are the same series with independent noise).
- **The GP samples**: +0.037 +- 0.028 (3 of 5 networks higher, an unclear gain) and an unstable prediction (median skill -0.37, 2 of 5). This agrees with the real IRMA result (switch-on up, switch-off not).
- A single real series predicts the held-out series rather well in this setup (median 0.63) and the augmentations do not improve on it; the prediction is stable in all 5 networks only with 3 or more real series.

## What this says

- The augmentations of GREMA only restate the uncertainty of the series that is given; they cannot supply the information that a **new, independent series with a different initial state or perturbation** gives. This is the same reason why the earlier simulation from a fitted model made things worse (`input/IRMA/README.md`).
- When the data are not enough, the options that add information are: **new series** (the first additional one is the most valuable: +0.10 AUROC), **prior knowledge** (a short candidate list), and **a structure that needs fewer parameters** (modules, shared parameters; not tested). For planning an experiment: two or three different series are worth far more than many time points or replicates of the same one.

## Limits

- 5 networks, one generator (the model of GREMA, 10 genes, 20 points, 3% noise); the noise level of the augmentation (0.03) and the GP variant (v2 with noise, 10 samples; the IRMA analysis used 20) were fixed in advance and not varied; other noise levels, other numbers of copies or a generator that is not the model of GREMA were not tried.
- One held-out series per network; the standard errors of the paired gains are 0.006-0.03, so the small gains of the augmentations (0.02-0.04) are at the limit of what can be told from the noise, but they are far from the 0.10-0.14 of the real series.
