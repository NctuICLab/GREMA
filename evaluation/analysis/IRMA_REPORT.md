# IRMA: perturbation sets (T2) and number of generations (T3)

HFODE, after PR #18, IRMA switch-on (16 time points) and switch-off (21), one replicate each, `-t 5`, `IRMA_no_prior_knowledge.txt`.
20 non-self directed pairs, 8 true (random AUROC 0.5, AUPRC 0.4). Values are mean AUROC / AUPRC over the runs; `conf` = CONFIDENCE_LEVEL,
`vote` = VOTE_SCORE. The sd of the AUROC over runs is 0.00-0.09, so differences below about 0.05 are noise.
Baseline = `irma_fixed` of `evaluation/analysis/results/final_results.tgz` (`-g 500`, 10 runs per setting).
The final_results.txt of all runs (and the step-1 progress files of the T3 runs) are in `evaluation/analysis/results/irma_t2_t3_results.tgz`
(see `evaluation/analysis/results/README.md`).

## T2. Built-in perturbation sets (`-perturb 9 -noise N`, `-g 500`, 5 runs per setting)

| Setting | switch-on conf | switch-on vote | switch-off conf | switch-off vote |
|---|---|---|---|---|
| `-f 1` baseline (10 runs) | 0.705 / 0.596 | 0.781 / 0.754 | 0.683 / 0.572 | 0.605 / 0.544 |
| `-f 1`, noise 0.05 | 0.702 / 0.593 | 0.781 / 0.743 | 0.673 / 0.557 | 0.598 / 0.544 |
| `-f 1`, noise 0.1 | 0.702 / 0.593 | 0.770 / 0.728 | 0.632 / 0.524 | 0.631 / 0.579 |
| `-f 3` baseline (10 runs) | 0.779 / 0.732 | 0.762 / 0.667 | 0.708 / 0.660 | 0.609 / 0.574 |
| `-f 3`, noise 0.05 | 0.713 / 0.617 | 0.765 / 0.699 | 0.702 / 0.593 | 0.633 / 0.602 |
| `-f 3`, noise 0.1 | 0.717 / 0.621 | 0.752 / 0.725 | 0.665 / 0.545 | 0.622 / 0.624 |

The switch-off vote looked a little better with noise 0.1 after 5 runs, so 15 more runs were made (20 in all):

| switch-off, noise 0.1 | conf | vote |
|---|---|---|
| `-f 1`, baseline (10) | 0.683 / 0.572 | 0.605 / 0.544 |
| `-f 1`, 20 runs | 0.629 / 0.521 | 0.602 / 0.547 |
| `-f 3`, baseline (10) | 0.708 / 0.660 | 0.609 / 0.574 |
| `-f 3`, 20 runs | 0.663 / 0.545 | 0.605 / 0.600 |

- The solver now uses the perturbation sets (the fix of `handoff/analysis`), but they do not improve the ranking: every difference is within noise, and
  the apparent gain of the switch-off vote disappears with 20 runs. conf is lower with noise (switch-off 0.045-0.054, switch-on `-f 3` about 0.06) and the runs differ more (sd of the AUROC 0.05-0.09 vs 0.00-0.04).
- This agrees with the hypothesis that independent noise per time point adds no information about the network.

## T3. Number of generations (`-g` 500 / 1000 / 2000, 10 runs per setting)

conf, AUROC / AUPRC:

| | switch-on `-f 1` | switch-off `-f 1` | switch-on `-f 3` | switch-off `-f 3` |
|---|---|---|---|---|
| 500 (baseline) | 0.705 / 0.596 | 0.683 / 0.572 | 0.779 / 0.732 | 0.708 / 0.660 |
| 500 (rerun) | 0.705 / 0.596 | 0.685 / 0.573 | 0.770 / 0.709 | 0.708 / 0.650 |
| 1000 | 0.702 / 0.593 | 0.708 / 0.600 | 0.737 / 0.655 | 0.699 / 0.609 |
| 2000 | 0.702 / 0.593 | 0.708 / 0.600 | 0.773 / 0.714 | 0.708 / 0.600 |

vote, AUROC / AUPRC:

| | switch-on `-f 1` | switch-off `-f 1` | switch-on `-f 3` | switch-off `-f 3` |
|---|---|---|---|---|
| 500 | 0.781 / 0.754 | 0.605 / 0.544 | 0.762 / 0.667 | 0.609 / 0.574 |
| 1000 | 0.775 / 0.739 | 0.627 / 0.547 | 0.772 / 0.724 | 0.648 / 0.620 |
| 2000 | 0.772 / 0.737 | 0.657 / 0.551 | 0.776 / 0.730 | 0.671 / 0.692 |

TP / FP of the voted network are about the same for all `-g` (4.0 / 1.0-1.8 for `-f 1`, 4.0-4.9 / 0.7-1.6 for `-f 3`). The two `-g 500` batches
(the baseline of `irma_fixed` and the rerun that was needed for the progress files) differ by 0.01 or less in conf AUROC.

- The conf results are stable between 500 and 2000 generations (all differences within 0.05); several settings give the same network in all 10 runs
  (conf sd 0.000). The switch-off vote rises with `-g` by 0.05-0.06 (`-f 1` 0.605 -> 0.657, `-f 3` 0.609 -> 0.671), at the border of the noise, and vote is the worse output on switch-off.

Convergence of the GA (`evaluation/analysis/scripts/convergence.py`, step-1 progress files, 5 genes x 30 GA runs x 10 runs): relative fitness gain
from generation 0.5g to 0.9g, median / 90th percentile:

| | `-g 500` | `-g 1000` | `-g 2000` |
|---|---|---|---|
| switch-on `-f 1` | 1.51% / 11.2% | 2.04% / 10.8% | 1.72% / 6.8% |
| switch-off `-f 1` | 0.61% / 4.0% | 0.80% / 4.2% | 0.24% / 1.6% |
| switch-on `-f 3` | 0.65% / 8.5% | 1.25% / 8.1% | 0.54% / 3.5% |
| switch-off `-f 3` | 0.48% / 3.5% | 0.50% / 3.3% | 0.14% / 1.4% |

- The fitness gain in the second half of the run is small (median below 2%); switch-on with `-f 1` converges slowest, but its network does not change.
- The share of GA runs that still improved in the last 0.2g generations is 84-91% for `-g 500` and 99-100% for 1000 / 2000: it depends on the window and says little; the gain is the better measure.
- 500 generations are enough on IRMA. The same check on SOS (`-g 500` / `-g 1000`) is in `evaluation/analysis/SOS_REPORT.md`.

## Correction to the hand-off text

`-f 3` is not clearly better than `-f 1` on IRMA: switch-on mean(conf, vote) AUROC is 0.778 vs 0.782 (`-f 1` vs `-f 3`), switch-off 0.618 vs 0.641, both within noise.
On SOS `-f 1` is clearly better, so `-f 1` is the safer default.

## Limits

- 5-20 runs per setting and 8 true edges; two IRMA series only. The T3 runs of `-g 500` (rerun) were made with the same code and options as the baseline, in a different session.
