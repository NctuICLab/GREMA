# When the measured variable is not the regulator: promoter activity, protein regulators, a shared driver (why SOS fails)

Question: on the SOS data the signs found by GREMA are wrong (every found LexA edge is `+`, LexA is a repressor) and the detection is weak, while dynGENIE3 scores higher. What does that depend on? Simulated networks with known edges and signs
(`evaluation/analysis/scripts/sim_protein_net.py`; archive `evaluation/analysis/results/sim_protein_results.tgz`; `python3 evaluation/analysis/scripts/eval_protein_sims.py <unpacked archive>` prints the tables). Earlier: signs of GREMA are reliable when the
measured variable is the regulator (`evaluation/analysis/SIGN_SIMULATION_REPORT.md`); a single-hub simulation with a hidden LexA protein reproduced the wrong signs of SOS (`evaluation/analysis/SOS_REPORT.md`).

## Simulation (rules fixed before the runs)

- 5 random signed networks (seeds 1-5; 10 genes, 1-2 regulators per gene, each edge `+` / `-` with probability 1/2; 14-18 edges, 80 in all), 4 series x 20 points, 3% noise, one global scaling to [0.01, 1].
- The regulators are **proteins** `P_j`; the **measured** variable is the **promoter activity** `a_i = TransMax_i prod_j 1/(1 + (k_ij/P_j)^(s_ij n_ij)) + b_i` of the regulated gene; `dP_j/dt = d_j (a_j - P_j)`, `d_j` in U(0.1, 0.3).
- Conditions: no driver (the series differ by random initial proteins); **transcription driver** (every promoter activity multiplied by `1 + u_s(t)`, `u_s(t) = s exp(-t/6)`, `s` = 0.5, 1, 1.5, 2 for the 4 series; all series start from the equilibrium; the proteins still follow the activity);
  **cleavage driver** (every protein is also cleaved, `-0.4 u_s(t) P_j`: the protein falls, the promoters it represses rise, as for LexA after UV).
- GREMA: default settings (3 slots, `-f 1`, `-g 500`), 1 run per network. dynGENIE3: RF and ET, 5 seeds.
- Collinearity of the measured genes (mean |correlation|) and how often the sign of the correlation of a true pair is its true sign: no driver 0.39 / 71%, transcription 0.56 / 54%, cleavage 0.64 / 69%.

## Results

**1. A shared driver and the opposite direction** (GREMA; random AUROC 0.5, AUPRC 0.18; the measured signal is the regulator):

| Condition | AUROC / AUPRC | true edges called as a regulation | of those, sign right | K best by vote: true edges |
|---|---|---|---|---|
| no driver | 0.770 / 0.443 | 63 / 80 (79%) | 62 / 63 (**98%**) | 32 / 80 |
| transcription driver (collinear, same direction) | 0.661 / 0.330 | 48 / 80 (60%) | 47 / 48 (**98%**) | 30 / 80 |
| cleavage driver (collinear, **opposite direction**) | **0.578** / 0.276 | 36 / 80 (45%) | 30 / 36 (**83%**) | 22 / 80 |

(When the measured variable is the regulator and the model is right: 0.914 / 0.711, signs 100%; `SIGN_SIMULATION_REPORT.md`.)
Collinearity lowers the detection (0.77 -> 0.66) and leaves the signs right; only when the protein moves opposite to the measured signal (cleavage) do the signs go wrong (83%) and the detection fall to random level.

**2. dynGENIE3 on the same data** (no signs; AUROC / AUPRC, 25 runs each): no driver RF 0.723 / 0.512, ET 0.746 / 0.562; transcription 0.646 / 0.383, 0.651 / 0.386; cleavage 0.609 / 0.306, 0.605 / 0.337.
It loses as much as GREMA (cleavage: 0.61 vs GREMA 0.58, within the noise): no sign that it copes better with collinearity or the opposite direction.

**3. Filtered regulator (negative result).** Idea: the regulator is the accumulated protein, so give GREMA a filtered copy of the measured signal as the only regulator (`augment_filtered.py`, `F(t+1) = F(t) + D (x(t) - F(t))`, `F(0) = x(0)`; the targets stay the measured genes;
`D` = 0.1, 0.2, 0.3 fixed beforehand; 3 slots; the pairs `F_j -> G_i` are read as `G_j -> G_i`, `unfilter_results.py`; no change of the solver). No driver data:

| Regulator | AUROC / AUPRC | called | sign right | K best: true edges |
|---|---|---|---|---|
| the measured signal | **0.770** / 0.443 | 63 / 80 | 62 / 63 (98%) | 32 / 80 |
| filtered, D = 0.1 | 0.688 / 0.388 | 41 / 80 | 39 / 41 (95%) | 28 / 80 |
| filtered, D = 0.2 | 0.727 / 0.481 | 43 / 80 | 41 / 43 (95%) | 34 / 80 |
| filtered, D = 0.3 | 0.742 / 0.408 | 50 / 80 | 46 / 50 (92%) | 36 / 80 |

No gain (the differences are below 0.05 or have different signs for AUROC and AUPRC): the measured signal already contains most of the information; a filtered regulator is **not recommended as a change of the solver**.

## What this says about SOS

- The wrong signs of SOS are best explained by the cleavage case: the measured promoter activity of lexA rises together with its targets while the protein falls (the correlation of the levels is 0.6-0.97), so GREMA, which reads the sign from the measured regulator, gives `+`.
  The simulation reproduces only part of it (83% right, not all wrong); the single-hub simulation (`SOS_REPORT.md`) gives all `+`; the strength of the driver (0.4, tau = 6) was not scanned; the real SOS data may have other causes (it is a single hub, 8 genes).
- dynGENIE3 does not solve this problem. Its high score on SOS does not come from the regulation: the experiment-permutation control leaves it unchanged (0.856 -> 0.822), the time index alone predicts dx/dt as well as all genes (R^2 0.49 vs 0.51), the best edges come from the "clock" genes G4, G2, G5,
  and the gold standard is a star around LexA (`METHODS_REPORT.md` section 1). It gives no signs, so the sign cannot be wrong; the regulator-level score of GREMA (0.914) is higher than that of dynGENIE3 (0.871).
- The way to get the regulation from such data is a model that knows the protein is cleaved (a protein state; a structural change, not tested) or data with perturbations (like Hackett), not a different edge scoring.

## Limits

- 5 networks, 1 GREMA run per network and condition, 80 edges; the sign accuracy of the cleavage condition rests on 36 called edges (6 wrong). The strength of the drivers and the decay rates were chosen once, no scan.
- The simulations are simple (HFODE-form promoters, a protein that is a smoothed activity); no feedback loops by design, no hub.
- Not tested: a model with a protein state, the filtered regulator on the driver data, SOS-like single hubs with random parameters.
