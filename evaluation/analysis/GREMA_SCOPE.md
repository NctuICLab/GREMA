# What GREMA is for, and what the evidence says

A short statement of scope, written from the analyses in this folder. Every line points to the report that has the numbers; reports of the open pull requests are named with their number (#30-#35) until they are merged.

## What GREMA does

GREMA builds a gene regulatory network from **time-series expression data**: for every gene it fits an ODE (a Hill-function model, `HFODE`, or an S-system) with a genetic algorithm and reports which regulators the gene needs, with a sign (the sign of the Hill coefficient), a confidence level and a vote score.
It takes **prior knowledge** (`+ - ? 0` for every pair), real **sampling times** (`-time`), an **external input** such as a stimulus (`-input`, #33) and a **data check** (`check_data.py`) that says before the run whether the data can support an inference. The result is a model that can be simulated.

## Use it when

- the data are a **perturbation time course** of a **small network** (5-60 genes) with **several series** that differ (different perturbations or initial states), **at least 8 time points** (the rule of thumb of `check_data.py`), and the measured variable is the **activity of the regulator** (mRNA of a TF that acts as mRNA);
- you can give a **candidate list or prior knowledge** (fewer regulator slots per gene gave fewer false regulations, `SSYSTEM_SLOTS_REPORT.md`, `HACKETT_REPORT.md`), and you want a **mechanistic, simulable model**, not only a ranking of pairs.

What it did on such data:
- **Simulations with the model correct**: AUROC 0.91 and signs 100% right with 4 series of 20 points; 0.82 and 92% with 1 series; 5 time points per series still 0.87 and 100% (`SIGN_SIMULATION_REPORT.md`).
- **Not only for Hill-shaped regulation**: with 4 series the fit tolerates an OR gate, non-monotone and linear regulation about as well as the Hill form (signs 93-100% right); a protein stage between the mRNA and the regulation lowers the detection by 0.13; with 1 series a wrong form costs 0.09-0.22 (`NONHILL_FORM_REPORT.md`, #34).
- **IRMA** (5 genes, 1 series, feedback loops): detection AUROC 0.78 (switch-on) and 0.61 (switch-off) against 0.60-0.63 and 0.49-0.57 for dynGENIE3 (`METHODS_REPORT.md`, section 2).
- **An input term** recovers a known stimulus in simulations (input edges AUROC 0.95, signs 15/15, held-out condition predicted with skill 0.92; `INPUT_TERM_REPORT.md`, #33).

## Do not rely on it when

- **One stimulus drives all the genes** (SOS, a differentiation time course): the scores do not come from the regulation (the experiment-permutation control leaves them unchanged; `METHODS_REPORT.md`, section 1).
- **Few time points and one series**: the internal hMSC series (RNA-seq and ATAC-seq, 1 series of 5 points) gave no usable network; the only signal was the co-expression of cell-cycle genes (`HMSC_FEASIBILITY_REPORT.md`, #30).
- You need the **signs of the regulations from real data** (below), or a **prediction of a new condition** (below).

## What the evidence does not support (yet)

| Claim | What was found |
|---|---|
| The signs are reliable on real data | **No.** On real IRMA 25-50% of the signs of the found true edges are right (`SIGNS_REAL_DATA_REPORT.md`); with the **true structure told** GREMA is right for 43-50% (chance), the sign is the same in every GA run and is the sign of a linear regression on the same structure (14 of 14), and it flips between switch-on and switch-off for 4 of 7 edges: the cause is the data, not the model (`SIGN_DIAGNOSIS_REPORT.md`, #35). |
| The fitted model predicts a new condition | **No.** Fitted on IRMA switch-on it predicts switch-off with a free-running skill of -4.7 against keeping the first value; it fits its own condition well (0.95) (`IRMA_PREDICTION_REPORT.md`, #31). |
| GREMA detects better than model-free methods | **Not shown.** Better than dynGENIE3 on IRMA, equal or lower on the Hackett TF-induction data (the one data set where the comparison is valid: 0.67-0.72 against 0.72-0.79, and a lower AUPRC) (`METHODS_REPORT.md`). |
| The input term helps on real data | **No evidence.** No improvement on IRMA (forecasts) and SOS (held-out experiment; the reversed input is as good); the tests have little power (`INPUT_TERM_REPORT.md`, #33). |
| The Hill form is what limits GREMA | **Not supported**: the data amount and the unmodelled protein stage are (`NONHILL_FORM_REPORT.md`, #34). |

## How to get a result that can be trusted

1. Run `check_data.py` before the run, and the experiment-permutation control (`permute_replicates.py`) when all series come from one stimulus.
2. Plan the experiment for GREMA: **3-4 or more series that differ** (different perturbations, doses or start times, an input that changes within a series), **8-12 or more time points**, **replicates**, the **regulator's activity** measured (or the protein), a **small candidate list**, and **one condition held out** to test the prediction.
3. Report the detection (which regulators), not the signs, unless the signs are checked against a perturbation of the regulator or are the same in every series.
4. Judge the result against a simple baseline (a correlation, a model without regulation) and against its own shuffled-data control.

## What is planned and not done

A stability check of the signs across series (a sign that differs between conditions is not reported as a sign); the ATAC footprint (or another prior) as a soft prior in the choice of the regulator slots; a delay between an input and the mRNA; module-level networks for data with many genes and few points; a guide to the number of series and time points needed (a simulation). None of these has been tested.
