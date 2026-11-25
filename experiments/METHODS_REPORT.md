# GREMA against dynGENIE3, and attempts to improve GREMA (IRMA, SOS, Hackett)

Archive of the runs and inputs: `experiments/results/methods_comparison_results.tgz` (see `experiments/results/README.md`).
All numbers are AUROC / AUPRC of the ranking of the TF -> gene pairs against the gold standard, mean over the runs; `vote` = VOTE_SCORE, `conf` = CONFIDENCE_LEVEL of GREMA
(HFODE, `-f 1`, `-g 500`). The gold standard defines the evaluated pairs and the metrics only; no setting was chosen with it. The rules of every comparison were written
before its results were looked at (the exceptions are marked). Differences below about 0.05 AUROC are noise.

## 1. Data sets and what they can show

| | IRMA | SOS | Hackett (this work) |
|---|---|---|---|
| series | 1 (switch-on / switch-off, 16 / 21 points) | 4 experiments, 50 points | 10 experiments (one TF induced in each), 7 points (15 min grid) |
| variables | 5 | 8 | 56 (10 TFs may regulate, 46 genes) |
| evaluated pairs / true edges | 20 / 8 | 56 (non-self) / 7 | 550 / 18 (gold network 3 of Ma et al. 2014) |
| random AUROC / AUPRC | 0.5 / 0.4 | 0.5 / 0.125 | 0.5 / 0.033 |

Hackett data: GEO GSE142864 (a **microarray**, not RNA-seq), design and rules in `experiments/HACKETT_REPORT.md`.

**SOS cannot rank methods by their regulation inference (control experiments).**
- A random forest predicts the rate of change dx/dt of the genes from the time index alone as well as from all other genes (out-of-bag R^2 0.49 vs 0.51; one gene alone: G1 0.41, G4 0.37, G5 0.34, G2 0.27, G6-G8 about 0): all genes follow the same UV time course, so any gene works as a "clock".
- Experiment permutation (`experiments/scripts/permute_replicates.py`: every gene takes its series from a randomly chosen experiment, time points stay aligned; this destroys the coupling of the genes inside an experiment and keeps the shared time course):

| | original | permuted |
|---|---|---|
| SOS, dynGENIE3 RF / ET (10 permutations) | 0.856 / 0.826 | 0.822 / 0.801 |
| SOS, GREMA vote (5 runs) | 0.750 | 0.720 |
| Hackett, dynGENIE3 RF / ET (10 permutations) | 0.723 / 0.791 | 0.548 / 0.557 |
| Hackett, GREMA vote, 2 slots (3 runs) | 0.736 | 0.516 |
| Hackett, Hill dictionary score (10 permutations) | 0.765 | 0.521 (sd 0.09) |

  On SOS the scores do not fall; on Hackett every method falls to random. So the Hackett scores come from the regulation inside the experiments, the SOS scores do not.
- The gold standard of SOS is a star around LexA (G2). Ranking the pairs by the mean score of the regulator alone (target ignored) gives 0.871 (dynGENIE3 RF), 0.914 (GREMA vote), higher than the edge scores (0.856, 0.751).
- Earlier statements that dynGENIE3 is better than GREMA on SOS therefore say little about regulation inference. (IRMA has one series per switch, the permutation cannot be made.)

## 2. dynGENIE3 against GREMA

dynGENIE3 (Huynh-Thu and Geurts 2018, original Python code, default parameters, RF and ET, 10 seeds, `experiments/scripts/dyngenie3_runner.py`; same input files as GREMA, the candidate regulators restricted like the knowledge file).

| Data | dynGENIE3 RF | dynGENIE3 ET | GREMA vote (conf) |
|---|---|---|---|
| IRMA switch-on | 0.599 / 0.579 | 0.631 / 0.579 | 0.781 / 0.754 (0.705 / 0.596) |
| IRMA switch-off | 0.492 / 0.414 | 0.572 / 0.516 | 0.605 / 0.544 (0.683 / 0.572) |
| SOS (not valid, see 1) | 0.856 / 0.393 | 0.826 / 0.313 | 0.751 / 0.317 |
| Hackett, GREMA 10 / 3 / 2 slots | 0.723 / 0.178 | 0.791 / 0.250 | 0.671 / 0.701 / 0.720 (AUPRC 0.058 / 0.080 / 0.101) |

GREMA is better on IRMA (one series, feedback loops, 15-20 time steps); dynGENIE3 is equal (RF) or better (ET, +0.07) on Hackett, and its AUPRC there is much higher (0.18-0.25 vs 0.06-0.11).
dynGENIE3 gives no signs. No evidence was found that GREMA beats dynGENIE3 on the data set where the comparison is valid.

## 3. What the edge scores of GREMA are

- With `-min 0` a candidate counts as used when it is in the random subset of regulator slots of a GA run (the Hill coefficient is exactly 0 otherwise; on SOS 56.6% of the entries are 0, 99% of the others have |n| > 0.1 and the median |n| is 3.97 with the bound 4): VOTE_SCORE is the fitness-weighted frequency of being in the subset.
- Contrast score (`experiments/scripts/contrast_scores.py`, fixed beforehand: Welch t of log fitness without vs with the candidate in the first-iteration runs): **worse** than vote: SOS 0.424 (vote 0.750), IRMA on 0.392 (0.785), IRMA off 0.713 (0.609), Hackett 3 / 2 slots 0.480 / 0.422 (0.707 / 0.736). Negative result.
- Post-processing of vote (`experiments/scripts/postprocess_scores.py`, fixed beforehand; paired runs): A1 = z score of vote within every target, A2 = vote / sum of the target, B = mean rank of vote and dynGENIE3 (RF or ET):

| Data | vote | A1 | A2 | dynRF | B (+RF) | B (+ET) |
|---|---|---|---|---|---|---|
| IRMA on | 0.785 | **0.829** | 0.785 | 0.599 | 0.661 | 0.727 |
| IRMA off | 0.609 | **0.706** | 0.618 | 0.492 | 0.582 | 0.621 |
| Hackett 3 slots | 0.707 | 0.701 | 0.695 | 0.729 | 0.754 | 0.776 |
| Hackett 2 slots | 0.736 | 0.767 | 0.736 | 0.729 | 0.767 | 0.785 |

  A1 helps on IRMA (+0.04 / +0.10) and slightly on Hackett with 2 slots (+0.03), nothing on Hackett with 3 slots; B stays between its two parts. No post-processing helps on all data.
- Hill coefficient bound and threshold (Hackett, 2 slots, 3 runs each): `-nmax 8` 0.726 / 0.115, `-min 0.5` 0.724 / 0.107, `-min 1` 0.728 / 0.108, against 0.736 / 0.107: nothing changes. The coefficients grow to the bound (median |n| 3.45 with the bound 4, 6.33 with 8), so their size carries no information and a threshold does not remove regulations (about 111 of the 3136 pairs are called with every setting). Negative result.

## 4. Candidate filters before GREMA

Three filters choose, for every target gene, the K = ceil(candidates / 2) best candidates (the others get `0` in the knowledge file), with 200 (LASSO), 100 (HSIC-Lasso) and 200 (Hill dictionary) stability-selection subsamples of 50% of the rows
(`experiments/scripts/lasso_prefilter.py`, `hsic_prefilter.py`, `hill_dictionary.py`; rules in the docstrings). Final score of a pair = vote of the GREMA run (0 for the filtered pairs) + 1e-4 x the filter score (`eval_prefilter.py`).

| Filter | model of the target's rate of change |
|---|---|
| LASSO | linear, sparse (own x_i projected out) |
| HSIC-Lasso | non-linear (kernel dependence), sparse |
| Hill dictionary | additive sum of Hill terms (18 per candidate: n in +-1, +-2, +-4, k at the 25 / 50 / 75 percentile), sparse; **additive, HFODE multiplies** |

Filter alone, and true edges kept (AUROC / AUPRC):

| Data | LASSO | HSIC-Lasso | Hill dictionary (3 seeds) |
|---|---|---|---|
| IRMA on | 0.750 / 0.714 (6 of 8 kept) | 0.646 / 0.531 (4) | 0.552 / 0.512 (5) |
| IRMA off | 0.792 / 0.661 (6) | 0.812 / 0.662 (6) | 0.594 / 0.540 (5) |
| Hackett | 0.642 / 0.056 (13 of 18) | 0.643 / 0.130 (12) | **0.765 / 0.132** (15) |
| SOS (not valid) | 0.176 / 0.084 (1 of 7) | 0.643 / 0.251 (4) | 0.398 / 0.128 (2) |

End to end (GREMA on the filtered knowledge, final score as above; unfiltered GREMA vote in the last column):

| Data | LASSO + GREMA | HSIC-Lasso + GREMA | Hill dictionary + GREMA | GREMA alone |
|---|---|---|---|---|
| IRMA on | 0.749 / 0.714 | 0.608 / 0.506 | 0.592 / 0.566 | 0.785 / 0.765 |
| IRMA off | **0.842 / 0.745** | 0.783 / 0.632 | 0.669 / 0.649 | 0.609 / 0.528 |
| Hackett, 2 slots | 0.697 / 0.093 | 0.664 / 0.171 | **0.768 / 0.127** | 0.736 / 0.107 |

(Hackett with a slot for every kept LASSO candidate, about 5 per gene: 0.673 / 0.073; SOS with HSIC-Lasso: 0.669 / 0.288 against 0.750 / 0.324.)

- The filter alone is often as good as the filtered GREMA (Hackett, Hill dictionary: filter 0.793 for seed 0, GREMA vote alone 0.741, combined 0.768): the genetic fit does not add to the filter there.
- IRMA switch-off improves with two filters (+0.17 / +0.23) and gets worse or stays with the others, so there is no filter that works everywhere; a filter that drops the true regulators (coverage) limits GREMA.
- On the first TF subset the deterministic Hill dictionary is the best filter (0.765, close to dynGENIE3 ET 0.791) and needs minutes instead of about one hour; as an additive model it is a replacement of the search, not an improvement of HFODE. **This was not reproduced on a second TF subset (0.698, section 7): the first subset was the favourable one.**

## 5. Why the data sets differ so much (partly tested)

- Samples against parameters: the Hill dictionary has 72 columns for about 15 samples on IRMA and 180 columns for about 60 samples on Hackett; on IRMA it is at random level (tested by the numbers).
- Type of signal: on Hackett the information is "gene i responds when TF j is induced inside the experiment" (checked by the permutation); on IRMA it is a time order within one series with feedback loops (not tested).
- Evaluation noise: 8 true edges on IRMA, a single edge moves the AUROC by more than 0.05; also switch-on and switch-off of the same network differ a lot (GREMA 0.785 vs 0.609).

## 6. Code and documentation changes in this PR

- `GREMA_main.pl`: the GA results of every iteration are kept (`<output>/iga_results/<gene>_step<k>.txt`, switch off with `GREMA_KEEP_ITER=0`); `<gene>.txt` was overwritten by each iteration. New option `-nmax` (bound of the Hill coefficients, default 4 as before).
- `README.md`: both are documented. New scripts in `experiments/scripts/` (see the docstrings): `dyngenie3_runner.py`, `permute_replicates.py`, `contrast_scores.py`, `postprocess_scores.py`, `lasso_prefilter.py`, `hsic_prefilter.py`, `hill_dictionary.py`, `eval_prefilter.py`.
- The dynGENIE3 code (https://github.com/vahuynh/dynGENIE3) was not copied into the repository; the Python environment was built with pixi (numpy, scikit-learn, scipy, pandas, pyHSICLasso).

## 7. Independent check with a second TF subset of the Hackett data

Same procedure (`hackett_prep.py --per-tf 10 --network 3 --step 15`), but the TFs ranked 11-20 by the number of responding genes (`--skip-tfs 10`; fixed before the runs): GCN4, SWI6, ROX1, RME1, SOK2, STB5, RTG3, NRG1, SKO1, CIN5
(55 variables, 540 pairs, 18 gold edges; the first subset was SFP1, GLN3, STE12, YOX1, PHO4, MOT3, YAP1, SUM1, ACE2, RLM1; the TFs are disjoint, some genes are shared). Archive: `experiments/results/hackett_set2_results.tgz`.
Compared (fixed beforehand): dynGENIE3 RF / ET (10 seeds), the Hill dictionary alone (3 seeds), GREMA with 2 slots (3 runs), the Hill dictionary filter + GREMA (3 runs), and the experiment-permutation control for dynGENIE3 and the Hill dictionary.

| Method | first TF subset | second TF subset |
|---|---|---|
| dynGENIE3 ET | **0.791** / 0.250 | 0.707 / 0.103 |
| dynGENIE3 RF | 0.723 / 0.178 | 0.701 / 0.075 |
| Hill dictionary alone | 0.765 / 0.132 | 0.698 / 0.083 |
| Hill dictionary filter + GREMA | 0.768 / 0.127 | 0.667 / 0.075 |
| GREMA, 2 slots (vote) | 0.736 / 0.107 | **0.609** / 0.050 (conf 0.602) |
| permutation control: Hill dictionary / dynGENIE3 RF / ET | 0.521 / 0.548 / 0.557 | 0.510 / 0.509 / 0.491 |

- The permutation control falls to random level on the second subset too: the scores come from the regulation inside the experiments.
- The three fast methods are equal within the noise on the second subset (0.70-0.71). The advantage of dynGENIE3 ET (+0.055) and the Hill dictionary over GREMA on the first subset is smaller or gone: only GREMA with the genetic search fell clearly (0.736 -> 0.609, the signs of its true edges are mixed: 6 `+`, 15 `-`).
- The Hill dictionary filter kept 13 of the 18 true edges; the filtered GREMA (0.667) is better than GREMA alone (0.609) and worse than the filter alone (0.703), as on the first subset: the genetic fit does not add to the filter.
- What can be said: on the Hackett data the fast methods (dynGENIE3, the Hill dictionary) are at 0.70-0.79, the genetic search of GREMA at 0.61-0.74, the ordering between the fast methods changes with the TF subset, and the AUROC of GREMA changes by 0.13. Two subsets of 18 edges, 3 GREMA runs each: the numbers are rough. Both subsets come from the same experiments.

## 8. Limits and corrections

- One valid data set with a larger gold standard (Hackett; 18 true edges, 3 GREMA runs per setting, gold network 3 is inclusive), one more (IRMA, 8 edges); SOS is not valid for comparing methods. The permutation control exists for SOS and Hackett only.
- Many settings were compared; the rules were fixed beforehand, all results (also the negative ones) are reported, but the multiplicity is not corrected.
- The Hill dictionary and the other filters were chosen on the first TF subset; the second subset (section 7) did not confirm the advantage seen there.
- Correction made during the work: a check of the Hill coefficients first used a column index shifted by one; the corrected numbers (56.6% of the entries are exactly 0, median |n| 3.97) give the same conclusion.
- The GREMA runs of the permutation control and the end-to-end tests are 3-5 runs per setting.
