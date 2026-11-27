# Is the limit of GREMA the Hill form or the data? Simulations with regulation functions that are not Hill products

Question: on real data the signs and the predictions of GREMA are poor (`SIGNS_REAL_DATA_REPORT.md`, `IRMA_PREDICTION_REPORT.md`). Is it because the true regulation is not of the Hill-product form that GREMA fits, or because of the data?
Archive of inputs and results: `evaluation/analysis/results/sim_nonhill_results.tgz` (see `evaluation/analysis/results/README.md`). Scripts: `sim_nonhill.py`, `eval_nonhill.py`. HFODE, `-f 1`, `-g 500`, no prior knowledge, default slots, 4 threads.

**Answer: in these simulations the Hill fit tolerates three of the four non-Hill forms when there are enough series; the form that hurts is the protein stage (a delay between the mRNA and the regulation), and few data amplify the damage of the others. The Hill form is not the main limit; the data and the unmodelled protein stage are.**
This is a limited test (one 5-gene topology, 6 seeds, mild non-Hill forms): it does not show that the Hill form describes real regulation.

## Design (the rules were fixed before the runs)

- Topology of IRMA (`sim_gen.py`): 8 regulations on 5 genes. Same seeds give the same parameters and initial values in every form, so the comparison is paired. 16 points, 5% noise, one global scaling to [0.01, 1], random initial values, a held-out series of the same system.
- Forms: **hill** (the product of Hill factors: the model of GREMA, the reference), **or** (the factors of the regulators are combined by an OR gate `1 - prod(1 - f)` instead of the product), **bell** (G2 -> G3 and G3 -> G5 are non-monotone, `4r/(1+r)^2` with `r = (x/k)^n`: the regulator activates at low and represses at high levels),
  **delay** (the regulators act through a protein that follows the mRNA with a rate 0.1-0.3, only the mRNA is measured), **linear** (production `clip(c + sum +- w x_j, 0, 2)`, no saturation).
- Data amount: **1 training series** (as in IRMA) and **4 training series**; 6 seeds each, 60 GREMA runs; dynGENIE3 (random forest, 5 seeds) on the same data as a model-free reference.
- Rules: AUROC drop of a form = AUROC(hill) - AUROC(form) with the same number of series. Drop < 0.05 and sign accuracy >= 85%: robust; drop >= 0.10 or sign accuracy < 70%: the form limits the Hill fit; otherwise moderate. The effect of the data is compared with the largest effect of a form.

## Results (mean +- sd over the 6 seeds; AUROC of the vote score against the 8 true regulations of the 20 pairs)

| Form | 4 series: AUROC GREMA | dynGENIE3 | sign right | 1 series: AUROC GREMA | dynGENIE3 | sign right |
|---|---|---|---|---|---|---|
| hill (reference) | 0.739 +- 0.081 | 0.773 | 25/26 | 0.701 +- 0.120 | 0.602 | 23/23 |
| or | 0.818 +- 0.046 | 0.748 | 27/29 | 0.613 +- 0.123 | 0.693 | 13/19 |
| bell | 0.808 +- 0.072 | 0.792 | 17/17 | 0.715 +- 0.077 | 0.689 | 18/18 |
| delay | **0.608 +- 0.102** | 0.710 | 18/18 | **0.484 +- 0.123** | 0.523 | 10/12 |
| linear | 0.725 +- 0.063 | 0.620 | 26/26 | **0.523 +- 0.121** | 0.516 | 10/14 |

Paired by seed (AUROC of the form minus AUROC of hill, mean +- standard error):

| Form | 4 series | 1 series |
|---|---|---|
| or | +0.079 +- 0.035 | -0.089 +- 0.057 |
| bell | +0.069 +- 0.032 | +0.014 +- 0.071 |
| delay | **-0.131 +- 0.059** (5 of 6 seeds lower) | **-0.218 +- 0.088** (5 of 6) |
| linear | -0.014 +- 0.039 | **-0.179 +- 0.036** (6 of 6) |

The effect of the data on the reference: 4 against 1 series +0.037 +- 0.041 (not clear).
By the rules: with 4 series or, bell and linear are *robust* and delay *limits* the Hill fit; with 1 series delay and linear (and or, by the sign accuracy 13/19) limit it.

- **Detection and signs with enough series**: the Hill fit loses nothing for the OR gate, the bell and the linear forms (or and bell are even higher than the reference, which is not a gain of the misspecification but a different information content of those data), and the signs of the called true edges stay right (93-100%).
  Only the protein stage lowers the detection (-0.13), and its signs are still right (18/18). The non-monotone edges of the bell form are found (vote 0.73 with 4 series, 0.43 with 1), their sign is not defined.
- **With one series the form matters**: linear (-0.18, 6 of 6 seeds), delay (-0.22) and, less clearly, the OR gate lose detection and some signs (or 13/19, linear 10/14). A misspecified form needs more data than a correct one.
- **Prediction needs the data most**: with one series the free-running prediction of a held-out series of the same system diverges even for the Hill form (skill -450 +- 1000 against 0.87 with 4 series); the fit of one series is not stable outside the training data.
  With 4 series the Hill model predicts better than a no-regulation ODE only for the Hill data (0.87 against 0.48) and for the protein stage (0.57 against -1.78); for or and linear the no-regulation ODE is as good (0.99 / 0.98, 0.84 / 0.80), for bell the Hill model is not better (0.55 against 0.56).
- **dynGENIE3** (no form assumed) is not clearly better or worse than GREMA on any form with 4 series (differences of 0.02-0.10 with sd 0.06-0.13).

## What this says

- It does not support "the Hill form fits only a few networks" as the main explanation. With 4 series a monotone or even non-monotone regulation, an OR instead of an AND and a linear regulation are all found about as well as the true Hill form; the signs stay right.
  This agrees with the linear model being as wrong as GREMA on the real IRMA signs (`SIGNS_REAL_DATA_REPORT.md`): the wrong signs there are not a matter of the form.
- It supports two other explanations: **(1) the amount of data** (one series: a misspecified form loses 0.09-0.22 in AUROC and the prediction diverges) and **(2) the regulation through an unmeasured protein stage** (-0.13 / -0.22), which is the case of many real regulations measured by mRNA only.
- It does not explain the signs of real IRMA: the sign was right in all the simulated forms (also delay), but only 25-50% on IRMA. Something that these simulations do not contain (confounding, real noise, an input, regulation that is not in the sign structure of the 7 edges, the sign table itself) is still unexplained.

## Limits

- 6 seeds per condition, one topology of 5 genes with 8 regulations, noise 5%, 16 points; the standard errors of the differences are 0.03-0.09, so differences of less than about 0.07 are not resolved.
- The forms are mild (an OR gate changes 3 genes, the bell 2 edges, the delay is one rate); a stronger misspecification (all regulations non-monotone, additive effects of many regulators, bistable switches) was not tested.
- The reference itself is modest (AUROC 0.74 with 4 series): the vote of 20 pairs with default slots and no prior knowledge. A prior knowledge or a larger number of series would raise it and could change which forms are tolerated.
- The held-out prediction compares with a no-regulation ODE only; the baselines of a sparse linear or a model-free regression were not used.
