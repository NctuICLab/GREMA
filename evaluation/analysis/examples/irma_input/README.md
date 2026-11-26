# Example: IRMA switch-on and switch-off with the medium as an external input

The template for a data set with an input (`-input`, see the section "External inputs" of the README of the repository): two series, an input row `U1` that is 1 in the series with galactose (switch-on) and 0.01 in the series without it (switch-off), the check of the data, GREMA with and without the input, and a forecast.

```shell
cd EMA_HFODE && make && cd ..                                   # the solver, once
bash evaluation/analysis/examples/irma_input/run_example.sh out_irma_input 4    # about 20-30 minutes with 4 threads
```

| Step | Command (inside `run_example.sh`) | What it does |
|---|---|---|
| 1 | `make_joint.py input/IRMA joint16.txt 16` | the two series as the replicates of one profile, the first 16 time points of each |
| 2 | `make_input.py joint16.txt table.tsv U1 --out joint16_input.txt` | adds the rows of `U1` from `table.tsv` (one line per series: start, level, baseline; `1 0 1 0` = series 1, from time point 0 the level 1; `2 0 0 0` = series 2, level 0), scaled to 0.01-1 |
| 3 | `check_data.py joint16_input.txt --input U1` | the input differs between the series (2 distinct profiles): no warning for it; the warning "only 2 series" stays |
| 4 | `perl GREMA_main.pl -i joint16_input.txt -k know_in.txt -o out -g 500 -f 1 -t 4 -input U1` | `know_in.txt` is `IRMA_no_prior_knowledge.txt` plus `U1  G<i>  ?` for every gene |
| 5 | `forecast.py joint16_input.txt runs 12` | fit on the first 12 points, forecast the last 4 (free-running, the input of the series given), against the model without the input and two baselines |

Numbers found with this example (5 GA runs each; see `evaluation/analysis/INPUT_TERM_REPORT.md`): U1 -> G3 is called `-` in all 5 runs (vote 1.00), the other genes (almost) never;
the forecast of the points 13-16 has an RMSE of 0.081 with the input and 0.077 without it (no benefit; the last points are a plateau), and with a fit on the first 8 points and a forecast of the points 9-16 0.113 against 0.109 (`forecast.py joint16_input.txt runs 8` after the runs `in8_r*`, `no8_r*`).

To use the same steps for your own data: replace `make_joint.py` by your profile, write `table.tsv` for your stimulus (a dose per series, the start time, optionally `--shape decay --halftime H` and `--times` for real sampling times), and run the check before GREMA.
The input is only identifiable when it differs between the series or changes within them; in this example it is a label of the condition (constant within a series), which is the weakest case.
