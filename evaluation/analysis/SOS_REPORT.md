# SOS DNA repair: results of the hand-off analyses (T1 and follow-ups)

HFODE, `-g 500` unless stated, 5 runs per setting, 64 pairs (8 true: LexA = G2 -> all 8 genes), knowledge file with every pair `?`
including the diagonal (`input/SOS/SOS_all_unknown_knowledge.txt`), `-t 4`. Values are mean AUROC / AUPRC over the runs (random 0.5 / 0.125);
`conf` = CONFIDENCE_LEVEL, `vote` = VOTE_SCORE. The sd of the AUROC over runs is 0.01-0.07, so differences below about 0.05 are noise.
Gold standard: no sign is stored, so a true edge with the wrong sign still counts as TP.
The final_results.txt of all runs are in `evaluation/analysis/results/sos_results.tgz` (see `evaluation/analysis/results/README.md`).

## 1. Real data, 4 experiments as replicates (T1)

| Setting | TP / FP | binary | conf | vote | SEN / SPE / ACC |
|---|---|---|---|---|---|
| `-f 1` | 6.0 / 21.8 | 0.680 / 0.197 | 0.693 / 0.213 | 0.755 / 0.306 | 0.75 / 0.61 / 0.63 |
| `-f 3` | 3.4 / 23.0 | 0.507 / 0.131 | 0.507 / 0.131 | 0.639 / 0.217 | 0.42 / 0.59 / 0.57 |
| `-f 1`, `-g 1000` | 5.0 / 20.8 | 0.627 / 0.168 | 0.655 / 0.197 | 0.748 / 0.271 | |
| `-f 3`, `-g 1000` | 4.0 / 22.2 | 0.552 / 0.139 | 0.554 / 0.142 | 0.660 / 0.202 | |

Paper (Table 5): TP 7, TN 42, FP 14, FN 1 (SEN 0.88, SPE 0.75, ACC 0.77). Without the 8 self pairs (56 pairs, 7 true): `-f 1` conf 0.674 / 0.198,
vote 0.751 / 0.317; `-f 3` conf 0.544 / 0.143, vote 0.670 / 0.243.

- `-f 1` is clearly better than `-f 3` (about 0.1-0.17 AUROC); on IRMA the two are within noise, so `-f 1` is the safer default. No rule that picks
  `-f` without looking at the gold standard was found: the noise level of the data (second-difference estimate / range: IRMA 0.11-0.12, SOS 0.08) does not explain it.
- 500 -> 1000 generations changes nothing beyond noise; the fitness gain from generation 0.5g to 0.9g has median about 1% (90th percentile 3-5%)
  (`evaluation/analysis/scripts/convergence.py`, now for any `-g`).
- Per iteration (`-f 1`): the first iteration decides about 85% of the pairs; FP are already there after it.
- Signs: of the 30 true edges found with `-f 1` only 6 are `-` (5 x LexA -> LexA, 1 x LexA -> G1), the rest are `+` although LexA is a repressor
  (`-f 3`: all `+`, and LexA -> LexA is never found). A sign-aware TP would be about 1.2 per run instead of 6.0.
- Combining conf and vote (mean or mean rank, `evaluation/analysis/scripts/combine_scores.py`, rules fixed beforehand) gives about the same AUROC as the better of the two
  (0.748 / 0.337 for `-f 1`) and a smaller sd, but no gain.

## 2. Why the signs are wrong

Promoter activity of lexA (G2) is measured, not LexA protein. After UV, LexA is cleaved, so all repressed promoters, including lexA's own, go up together.
In the data the levels of G2 and the other genes are positively correlated (-0.12 to 1.0 depending on the experiment and gene, mostly 0.6-0.97), while G2 vs the derivative of the other genes is
negative (about -0.4 to -0.7).

- Fixing LexA -> LexA to `-` in the knowledge file (`-f 1`, 5 runs; a diagnostic that uses prior knowledge, not a setting): TP 4.8 / FP 21.0,
  vote 0.749 / 0.357, conf 0.645 / 0.197. The downstream edges stay `+`. So a missing self-repression term is not the cause.

## 3. Simulated SOS, hidden LexA protein (`evaluation/analysis/scripts/sim_sos.py`, `sim_sos_p.py`)

Hill repression of 8 promoters by a hidden LexA protein P (P represses its own promoter too), UV cleaves P (doses 5, 5, 20, 20), 4 replicates, 50 points every
6 min, 3% noise, scaled to [0.01, 1] like the real data. The output is the promoter activity, as in the real data. Seeds 1-5, `-f 1`, default settings.

| Data | pairs (true) | TP / FP | binary | conf | vote | signs of the TP edges |
|---|---|---|---|---|---|---|
| P hidden (`sim_sos.py`) | 64 (8) | 2.4 / 22.8 | 0.446 / 0.134 | 0.446 / 0.131 | 0.429 / 0.142 | LexA -> LexA `-` in 5/5 runs; the 7 downstream edges found are all `+` |
| P measured as G9 (`sim_sos_p.py`), `-f 1` | 72 (8: P -> G1..G8) | 2.0 / 19.0 | 0.477 / 0.133 | 0.465 / 0.130 | 0.658 / 0.246 | all 10 `-` |
| P measured, `-f 3` | 72 (8) | 0.6 / 22.6 | 0.361 / 0.108 | 0.353 / 0.108 | 0.535 / 0.165 | all 3 `-` |

- With P hidden the model also fails and the downstream edges it finds are `+`, as in the real data: a regulator that is not the measured variable is enough to explain
  the wrong sign. With P measured the sign is right (`-`).
- The accuracy is still low with P measured: all promoters are driven by the same signal and are strongly collinear (correlation about 0.97, stronger than in the real data),
  so edges between the promoters are found as well (about 135-146 per 5 runs, all false). So the signs are a data/measurement problem, but the number of edges found is not.
- `-f 3` is worse than `-f 1` also here (5 seeds, sd 0.04-0.15: the direction agrees with the real data, nothing more).
- 5 seeds only; the simulated curves are more similar to each other than the real ones (derivative correlation about -0.15 instead of -0.7).

## 4. Regulator slot cap (`GREMA_SLOT_DIV`, `EMA_HFODE/src/Getknowledge.cpp`)

With all regulators unknown, each chromosome holds only `round(possible / 3)` regulator slots (3 for SOS). The divisor is now read from the environment variable
`GREMA_SLOT_DIV` (default 3, so nothing changes without it). `-f 1`, `-g 500`, 5 runs:

| slots per gene | TP / FP | binary | conf | vote | SEN / SPE / ACC |
|---|---|---|---|---|---|
| 8 (`GREMA_SLOT_DIV=1`) | 6.4 / 33.6 | 0.600 / 0.155 | 0.630 / 0.175 | 0.704 / 0.239 | 0.80 / 0.40 / 0.45 |
| 3 (default, T1) | 6.0 / 21.8 | 0.680 / 0.197 | 0.693 / 0.213 | 0.755 / 0.306 | 0.75 / 0.61 / 0.63 |
| 2 (`GREMA_SLOT_DIV=4`) | 3.4 / 14.8 | 0.580 / 0.152 | 0.587 / 0.161 | 0.735 / 0.364 | 0.42 / 0.74 / 0.70 |

FP fall with fewer slots but so do TP; the default 3 is the best compromise for AUROC. The paper's FP (14) is reached at 2 slots only with TP 3.4.

## Limits

- 5 runs per setting, 8 true edges; one real data set. The simulated data are made by us and only show what is sufficient to produce the observed behaviour.
- How the paper fed the 4 experiments to GREMA and its preprocessing are not stated (we use `log2(x+1)`, one global linear scaling, 4 replicates), so part of the gap to the paper may come from that.
- The IRMA runs of the same hand-off (T2 perturbation, T3 `-g` 500/1000/2000) are not included here.
