# Signs and detection of GREMA on simulated signed networks, with the model correct and misspecified

Question: the signs of the regulations found on real data were unreliable (SOS: every edge `+` although LexA is a repressor; Hackett second TF subset: 6 `+` and 15 `-` among the found true edges). Is the sign mechanism of GREMA
(the sign of the Hill coefficient n, `+` for activation) itself reliable? Archive of the inputs and results: `evaluation/analysis/results/sign_simulation_results.tgz`
(`python3 evaluation/analysis/scripts/eval_sim_variants.py <unpacked archive>` gives the table).

## Simulation (`evaluation/analysis/scripts/sim_signed_net.py`; rules fixed before any run)

- Network: 10 genes, every gene has 1 or 2 regulators drawn from the other genes, every edge `+` or `-` with probability 1/2 (14-18 edges per network; over the 5 networks 34 of 80 edges activate, so "all +" would be right in 43%, random 50%).
- Dynamics in the form of the solver (Euler, step 1): `dx_i/dt = TransMax_i * prod_j 1/(1 + (k_ij/x_j)^(s_ij n_ij)) + b_i - degrade_i x_i`, `s_ij` = the sign, n in U(2, 4), k in U(0.2, 0.6), TransMax in U(0.6, 1.0), b in U(0, 0.05), degrade in U(0.2, 0.5).
  4 series of 20 points, random initial values (U(0.05, 0.9)), 3% multiplicative noise, one global scaling to [0.01, 1]. 5 networks (seeds 1-5; the same seed gives the same network in all variants).
- GREMA: default settings (3 slots per gene, `-f 1`, `-g 500`, all pairs `?`, diagonal `0`); baseline 2 runs per network, variants 1 run per network.
- Measures: AUROC / AUPRC of VOTE_SCORE; share of the true edges called as a regulation (`+` or `-`); among those the share with the true sign (random: 0.5; `eval_signs.py`, `eval_sim_variants.py`);
  the same for the K best pairs by vote (K = number of true edges).
- **The generating model of the baseline is the model of GREMA: this is an upper bound; real data are never so favourable.**

## Results (GREMA)

| Setting | runs | AUROC / AUPRC | true edges called | sign right | top-K by vote: true edges | sign right |
|---|---|---|---|---|---|---|
| baseline (the model is right) | 10 | 0.914 / 0.711 | 148 / 160 (92%) | **148 / 148 (100%)** | 109 / 160 (68%) | 109 / 109 |
| noise 10% | 5 | 0.879 / 0.628 | 72 / 80 (90%) | 71 / 72 (99%) | 50 / 80 (62%) | 50 / 50 |
| noise 20% | 5 | 0.855 / 0.585 | 68 / 80 (85%) | 67 / 68 (99%) | 46 / 80 (58%) | 45 / 46 |
| 8 time points per series | 5 | 0.903 / 0.673 | 74 / 80 (92%) | 74 / 74 (100%) | 51 / 80 (64%) | 51 / 51 |
| 5 time points per series | 5 | 0.873 / 0.643 | 70 / 80 (88%) | 70 / 70 (100%) | 50 / 80 (62%) | 50 / 50 |
| 1 series only | 5 | 0.815 / 0.517 | 65 / 80 (81%) | 60 / 65 (**92%**) | 45 / 80 (56%) | 43 / 45 |
| regulators act as the mean of their Hill terms (not the product) | 5 | 0.876 / 0.613 | 71 / 80 (89%) | 69 / 71 (97%) | 44 / 80 (55%) | 43 / 44 |
| 2 genes that are not measured (hidden) take part | 5 | 0.923 / 0.677 | 61 / 64 (95%) | 61 / 61 (100%) | 43 / 64 (67%) | 43 / 43 |

Random AUPRC 0.14-0.18. dynGENIE3 (no signs; RF / ET, 5 seeds, same data): AUROC 0.845 / 0.863, AUPRC 0.699 / 0.702 on the baseline networks, i.e. GREMA is better at detection here (0.914 / 0.711) because its model is the generating model.

## Conclusions

- **The sign mechanism of GREMA is reliable**: 148 of 148 called true edges have the right sign in the baseline, 99-100% with noise up to 20%, with 5 time points, with a different way of combining the regulators and with hidden genes; the lowest value (92%) is with a single series.
- **Detection suffers more than the sign**: AUROC falls from 0.914 to 0.82-0.88 and about one third of the true edges is not among the K best pairs. The number of series matters most (1 series: AUROC 0.815, 4 series: 0.914), as seen on IRMA and SOS.
- So the wrong signs on real data (SOS, Hackett) do not come from noise, few time points, the form of the combination or a plain hidden gene. More probable (not tested): the measured variable is not the regulator (SOS: promoter activity, the protein is cleaved) and strongly collinear genes.

## Limits

- Upper bound: the model of the baseline is the generating model; the misspecifications are mild (noise up to 20%, at least 5 points, 1-2 regulators per gene, a hidden gene that is an ordinary regulated gene, no strong collinearity, no feedback loops by design).
- 5 networks, 80 edges per variant, one run per network in the variants; the sign accuracy is counted among the true edges that were called a regulation (8-19% of the true edges were not).
- Not tested: a measured variable that is not the regulator (like the SOS promoter activity), feedback loops, larger networks, the sign on real data with a known sign.
