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
| Option | `-input U1,U2` names the rows that are inputs (not fitted; HFODE only; error with `-m s-system`, with `-perturb`/`-noise` they are copied unchanged into the extra series) |
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

1. **Simulation with a known truth** (`sim_irma_like.py` extended): HFODE data with a stimulus on two genes (one activation, one repression), 3 conditions with different levels or timing, 3 series of 16 points, 5% noise, 10 parameter sets. For each: gene-to-gene AUROC and sign accuracy with and without the input term, AUROC and sign accuracy of the input edges, free-running skill on a held-out condition (baselines: persistence, training mean, no-regulation ODE).
   Controls: a constant input; the input levels shuffled between the series.
   *Pass*: input-edge AUROC >= 0.8 and sign accuracy >= 90%; gene-to-gene AUROC not more than 0.05 below the fit without input; the held-out skill above the no-regulation ODE by at least 0.1; the two controls give input-edge AUROC near 0.5 or a warning of the data check.
2. **IRMA, switch-on / switch-off** (the test of PR #31 again): input = the medium (a step 0.01 to 1 for switch-on, 1 to 0.01 for switch-off), both conditions as two series (switch-off cut to the 16 points of switch-on, or the real times with `-time`); held-out cross-condition prediction both ways. *Pass*: free-running skill above the no-regulation ODE (0.44 on -> off, 0.64 off -> on) with the default slots.
3. **SOS, held-out dose**: input = the UV dose (4 experiments: 2 at 5 and 2 at 20 J/m^2, a step at the time of the exposure, 0.01 to a level proportional to the dose); fit on one dose, predict the other, both ways; the experiment-permutation control with the dose labels shuffled. The edges and signs (all LexA edges are `+` now) are secondary. Two input shapes (a step and a decay) to see how much the result depends on the shape.
   *Pass*: held-out skill above the baselines and the shuffled control clearly worse.

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
