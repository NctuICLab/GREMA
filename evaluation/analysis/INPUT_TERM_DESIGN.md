# Design: an external input for HFODE (`-input`)

Status: **design for review, nothing is implemented.** Decisions that need an answer are in section 7.

## 1. Why

GREMA fits `dx_i/dt` as a function of the other genes only. A stimulus (UV in SOS, the medium in IRMA, a ligand, a differentiation signal) cannot be written into the model, so (a) a stimulus that drives many genes at once is learned as regulation between the genes (SOS: the experiment-permutation control leaves the scores unchanged, every LexA edge is `+`), and (b) the fitted model does not carry over to a new input: fitted on IRMA switch-on it predicts switch-off with a free-running skill of -4.7 against persistence (0.44 for a model without regulation), `IRMA_PREDICTION_REPORT.md` (PR #31).
An input term would let the network explain only what the stimulus does not, give the stimulus-to-gene effects (which genes respond directly, with which sign), and allow the model to be simulated under a new input.

## 2. The model

For gene `i`, with regulating genes `R_i` and inputs `M_i` (both given by the knowledge file, as now):

```
dx_i/dt = TransMax_i * prod_{j in R_i} 1/(1+(k_ij/x_j)^n_ij) * prod_{m in M_i} 1/(1+(k_im/u_m(t))^n_im) + b_i - degrade_i * x_i
```

- An input is treated as a regulator whose trajectory `u_m(t)` is **given, not fitted**: the same Hill factor, the same sign convention (`n > 0` activation by the input, `n < 0` repression), the same parameter ranges, the same slots, the same vote score, the same `+ - ? 0` in the knowledge file.
  The existing fit already uses the observed values of every other variable (divide and conquer), so for the fit an input is simply a variable that is never a target.
- `u_m` is on the same scale as the data, in [0.01, 1] (GREMA cannot use 0): a step from 0.01 to 1, a dose as a level between 0.01 and 1, or a measured signal scaled with the same linear map as the genes. With `-time` an input is interpolated linearly inside an interval, like every other variable.
- Output: `final_results.txt` gets lines `U1  G3  +  confidence  vote` (input, gene, sign, ...), next to the gene-to-gene lines.

## 3. Interface

| Part | Change |
|---|---|
| Data file | an input is a row like a gene (`rep<k>  U1  v1 v2 ...`, same number of time points); a different input level in every replicate (series) is what makes it identifiable |
| Option | `-input U1,U2` names the rows that are inputs (not fitted; HFODE only; error with `-m s-system`, not with `-perturb` in the first version) |
| Knowledge file | `U1  G3  ?` etc. (input in the TF column); an input is never a target (a line with an input as GENE is an error); without a line an input does not regulate (`0`) |
| Variables | `G`s first, then inputs, so the gene numbers and `run.log` are unchanged; `num_var = num_tf = genes + inputs`; the genes that are fitted are the first `genes` (the loop of `GREMA_main.pl` runs over them only) |
| Slots | the candidates of a gene include the inputs marked `?` (rules and `GREMA_SLOT_DIV` as now) |
| Simulation | new `scripts/grema_simulate.py` (from `predict_irma.py`): the best row of each gene, or the ensemble of the runs, a data file with the initial values and the input series of a new condition -> predicted time courses |
| Input builder | new `scripts/make_input.py`: from a small table (one line per series: name, start time, level, duration, optional end level; or a column of measured values) it writes the input rows, scaled with the same linear map as the genes (`--like DATA.txt`) and aligned to the time points of the data (`-time`), ready to be appended to the data file; the user does not write the numbers by hand |
| Data check | `check_data.py` (with the data file, `-input` and the knowledge file) warns when an input is constant, when it has the same shape as the time index in every series (not identifiable), when it is identical in all series, when it is not on the scale of the genes, or when the series have different numbers of points; it also prints a design hint (at least two different levels or start times, at least 8 points) |
| Shape sensitivity | `scripts/input_sensitivity.sh` (or `make_input.py --shapes step,decay`): the same run for two shapes of the same input (a step and a step followed by a decay with a given half-time); a table of the input edges, their signs and the held-out skill for each shape, so that the user sees how much the result depends on the shape she chose |
| Example | `evaluation/analysis/examples/irma_input/`: IRMA switch-on and switch-off with the medium as input, from the original files to the held-out prediction, with the commands and the expected numbers; the template for new data |

Without `-input` nothing changes: the same input, config and knowledge files and the same solver command lines (the runs themselves are random, there is no seed option, so a regression test compares the files and the statistics of the results, not the bits).

Most of the work is in `GREMA_main.pl` (reading, ordering, loops, knowledge, output); the solver should need little or none, because variables that are only regulators already work like observed regulators (`num_tf` and the divide-and-conquer fit).
To be confirmed in the first step of the implementation: how `EMA_main.cpp` / `Model.cpp` read the data and the knowledge matrix when `num_var` is larger than the number of fitted genes, and the limit of 52 slots.

## 4. Alternatives

| Form | For | Against | Decision |
|---|---|---|---|
| **Input as a Hill regulator in the product (above)** | reuses everything (parameters, slots, sign, vote); saturating; small change | an input of 0.01 gates an activated gene to nearly 0 (the baseline `b` carries it); one more `k`, `n` per input and gene | **first version** |
| Additive input `+ beta_i u(t)` | one parameter; no gating | no saturation; a second kind of regulator in the code, in the chromosome and in the report | optional later, `-input-mode add` |
| Filtered input (a delay: `du_f/dt = (u - u_f)/tau`) | a protein stage between the signal and the mRNA, which lowers the detection in the simulations | a state and a parameter more | **phase 2**, after the first version has been tested |
| Fitted input shape | no need to know the profile | many free parameters, not identifiable with a few series | not planned |
| Time-varying parameters | changes of the network itself (chromatin) | out of scope | not planned |

## 5. Validation (criteria fixed before any run)

1. **Simulation with a known truth** (`sim_input_net.py`, `eval_input_sim.py`): HFODE data with a stimulus on two genes (one activation, one repression), 3 conditions with different levels or timing, 3 series of 16 points, 5% noise, 10 parameter sets. For each: gene-to-gene AUROC and sign accuracy with and without the input term, AUROC and sign accuracy of the input edges, free-running skill on a held-out condition (baselines: persistence, training mean, no-regulation ODE).
   Controls: a constant input; the input levels shuffled between the series.
   *Pass*: input-edge AUROC >= 0.8 and sign accuracy >= 90%; gene-to-gene AUROC not more than 0.05 below the fit without input; the held-out skill above the no-regulation ODE by at least 0.1; the two controls give input-edge AUROC near 0.5 or a warning of the data check.
2. **IRMA, switch-on / switch-off** (changed during the implementation): the held-out cross-condition prediction of PR #31 **cannot be repeated with an input**: in IRMA the medium is constant within each series (a label of the condition: galactose present in switch-on, absent in switch-off), so a model fitted on one condition sees only one value of the input and cannot learn its effect. Instead: both conditions as two series (switch-off cut to the 16 points of switch-on), input = the medium (1 in switch-on, 0.01 in switch-off); (a) the input edges and signs, and the fit of the two series with and without the input, (b) a forecast of the last 4 points of both series from a model fitted on the first 12, with and without the input, against the training mean and the no-regulation ODE. *Pass* (b): the free-running skill with the input above the model without the input and above the no-regulation ODE; (a) is descriptive. The held-out-condition test needs an input that changes inside a series, which is the case of SOS (test 3).
3. **SOS, held-out experiment** (made precise during the implementation, before any SOS run): input = the UV dose, raw 5 for experiments 1 and 2 and 20 for experiments 3 and 4, scaled as u = 0.01 + 0.99 raw / 20 (0.2575 and 1; `prep_sos_input.py`). The exposure is at the first time point, so a **step** is constant within every series (a label of the dose, like the medium in IRMA). Fitting on one dose and predicting the other (the first plan) has the same input in both training series: it has the shape of the time index and `check_data.py` calls it not identifiable. The **main design is therefore leave-one-experiment-out**: fit on three experiments (both doses present), predict the fourth, for each of the four experiments. The input is a **decay** (the damage is repaired): the dose at the start, falling to 0.01 with a half-time of **60 min (10 time points)**. Compared: the model with the input, the model without it (the experiments as replicates, as now), the training mean, the no-regulation ODE, and the control rev10 (the same decay in the reverse order of time: a stimulus that grows). 3 GA runs per fold and mode (36 runs). The single-dose design (`to20`, `to5`, decay and no input, 3 runs) is run as well and reported as non-identifiable. Sensitivity to the shape (half-times 30 and 120 min, the step; 2 runs per fold) is reported, not used to choose. The edges and signs (all LexA edges are `+` now) are secondary.
   *Pass*: the held-out skill (mean over the four folds) of the decay model above the model without the input and above the no-regulation ODE, and rev10 worse than the decay model by at least 0.1.
4. **The tools** (`make_input.py`, `check_data.py`, `input_sensitivity.sh`): `make_input.py` writes the same input rows as the hand-written ones of tests 1-3; `check_data.py` raises the warning for each of the four bad inputs (constant, same shape as the time, identical in all series, wrong scale) and none for the inputs of tests 1-3; the example in `examples/irma_input/` runs from a clean checkout with the commands of its README and gives the numbers in it.

Passing 1 only shows that the implementation works; the claim "the model carries over to a new input" needs 2 or 3 as well. A failure of 2 and 3 with a pass of 1 would mean the real systems differ from the model (the protein stage, for example), and is reported as such.

## 6. Limits that stay

**When an input does not help** (the data check says so): observational series without a stimulus; a stimulus that is a single step in one series (it has the shape of the time index, so it cannot be told from the program of the cell: this is the case of the present hMSC differentiation series); an input that is the same in every series. A useful experiment has at least two series with a different level or start time of the stimulus, and at least 8 time points.

The number of time points and series (`check_data.py`), signs on real data (25-50% on IRMA, cause unknown), regulation at the protein level when the mRNA is measured, large networks (module level and priors are separate work), and rules of the network that change with the cell state.

## 7. Decisions needed

1. **Form**: the input as a Hill regulator in the product (recommended, first version) or the additive form first?
2. **Input profile**: given by the user, scaled to [0.01, 1] (recommended; two shapes tested for SOS) or fitted?
3. **First scope**: HFODE only, with `-time` supported, S-system not (recommended, as for `-time`)?
4. **IRMA joint fit**: the two conditions have 16 and 21 points and the data format needs one number of points; cut switch-off to 16 points (recommended for the test) or allow series of different length (a larger change)?
5. **Tools for the user** (`make_input.py`, the extended `check_data.py`, `input_sensitivity.sh`, the example): in the first version, or the first version with only the format and the example?
6. **Phase 2** (a delay for the protein stage): in or out of this work?
