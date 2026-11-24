# SOS DNA repair system of *E. coli* in the GREMA format

The real data set used in the GREMA paper (Tsai et al., Bioinformatics 2020, Section 2.1 and 3.5): 8 genes
(uvrD, lexA, umuD, recA, uvrA, uvrY, ruvA, polB = G1..G8), promoter activities (GFP/min/OD, reporter plasmids)
after UV exposure, 50 time points with a 6-minute step.

**Source of the data.** Ronen et al. (2002), *Assigning numbers to the arrows: parameterizing a gene regulation
network by using accurate expression kinetics*, PNAS 99:10555. The original file ("Data on SOS system",
`sosdata.zip`: `Exp1.txt` ... `Exp4.txt`) is on the download page of Uri Alon's lab
(https://www.weizmann.ac.il/mcb/UriAlon/download/downloadable-data). It has 4 experiments: 1 and 2 with UV = 5 J/m^2,
3 and 4 with UV = 20 J/m^2. The Bioconductor package SPEM (`data/sos.rda`) holds only experiment 3, as `log2(x + 1)`.

**True network (gold standard).** The paper gives it only in Figure 5 (no list in the text or the supplement):
LexA represses all 8 genes, including itself (7 LexA -> other genes and the autoregulation LexA -> LexA, 8 regulations;
in the paper GREMA finds 7 of them and misses LexA -> uvrA). The paper counts the 64 gene pairs including the
self pairs (8 true, 56 false).

| File | Content |
|------|---------|
| `SOS_4exp_timeseries_expression.txt` | all 4 experiments as 4 replicates (`repeat_number=4`), 50 time points |
| `SOS_timeseries_expression.txt` | experiment 3 only (from SPEM), 1 replicate |
| `SOS_all_unknown_knowledge.txt` | knowledge file with every pair `?`, **including the diagonal** (as in the paper: 64 pairs) |
| `SOS_no_prior_knowledge.txt` | knowledge file with the diagonal fixed to `0` (56 pairs; `lexA -> lexA` cannot be found) |
| `SOS_GoldStandard.tsv` | gold standard, `TF GENE 0/1`, G2 (lexA) -> all 8 genes = 1; the sign (repression) is not stored |
| `convert_sos_raw.py` | original `Exp1..4.txt` -> `SOS_4exp_timeseries_expression.txt` |
| `export_sos.R`, `convert_sos.py` | SPEM `sos.rda` -> CSV -> `SOS_timeseries_expression.txt` and `SOS_no_prior_knowledge.txt` |

**Preprocessing.** `log2(x + 1)` (the transformation of SPEM) and then one linear map of all values to [0.01, 1]
(GREMA cannot use 0). The paper does not describe its preprocessing, and it does not say how its 4 experiments
were given to GREMA; treating them as 4 replicates is an assumption.

```shell
python3 input/SOS/convert_sos_raw.py <dir with Exp1.txt..Exp4.txt> input/SOS
mkdir -p output/SOS_4exp
perl GREMA_main.pl -i input/SOS/SOS_4exp_timeseries_expression.txt -k input/SOS/SOS_all_unknown_knowledge.txt \
    -o output/SOS_4exp -g 500 -f 1 -t 4
python3 experiments/scripts/eval_runs.py --gold input/SOS/SOS_GoldStandard.tsv --self-pairs \
    "run=output/SOS_4exp/final_results.txt"
```

For the self pairs to be evaluated (`--self-pairs`) the knowledge file must have `?` on the diagonal. The paper reports
TP 7, TN 42, FP 14, FN 1 over the 64 pairs (SEN 0.88, SPE 0.75, ACC 0.77) and per GREMA iteration (Table 5).
