# EMA_ssystem
This code is the Evolutionary Modeling Algorithm to solve the S-system model

## Notes on the search (IGA)
* Options: `-i gene -I iteration -n runs -G generations -P population -F objective -c cc_mode -s seed`
  (`-s` is optional; without it the seed is `time ^ (pid << 16)`, so processes started in the same second differ).
* `-G n` runs exactly `n` generations. `nochange k` (config file) stops a run early when the best fitness did not change
  for `k` consecutive generations; `k <= 0` switches the early stop off.
* Regulators of the target gene. Knowledge row = target gene, `?` free, `+` g>0, `-` g<0, `0` none.
  `max_reg = round(N/3)` limits only the additional `?` regulators: every `+`/`-` regulator is always in the model
  (the number of regulators is `max(max_reg, #confirmed)`; a row without `?` keeps all its regulators).
  The first slots of the chromosome are dedicated to the confirmed genes, so initialisation, crossover, mutation and
  decoding cannot lose them or give them the wrong sign. The other slots hold distinct `?` genes; location mutation
  and crossover repair never create duplicates.
* Output: `IGA_result<gene>_step<I>.txt` has one line per run: alpha, g(1..N), beta, h(1..N), delta_t, fitness, cc.
  `cc` is the correlation of observed and simulated profile of the best model of that run over all replicates
  (a constant profile gives 0). `calprofile<gene>_step<I>.txt` holds the profile of the best model of every run.
* Observed data have no point after the last time point, so the derivative term of the objective is skipped there.
