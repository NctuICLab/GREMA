# GREMA

**GREMA** (**G**ene networks **R**econstruction using **E**volutionary **M**odelling **A**lgorithm) is a program for inferring  a novel type of gene regulatory network (GRN) with confidence levels for every inferred regulation, which is emulated GRN (eGRN). The higher the confidence level, the more accurate the inferred regulation. GREMA gradually determines the regulations of an eGRN with confidence levels in descending order using either an S-system or a Hill function-based ordinary differential equation model. It makes use of an evolutionary modelling algorithm (EMA) that is based on evolutionary intelligence, including both crowd wisdom and an evolutionary strategy, to cope with the underdetermined problem. EMA uses an intelligent genetic algorithm to solve the large-scale parameter optimisation problem.

## Citing GREMA

GREMA: Modelling of emulated gene regulatory networks with confidence levels based on evolutionary intelligence to cope with the underdetermined problem. _Bioinformatics_, 2020, [PMID:32399550](https://pubmed.ncbi.nlm.nih.gov/32399550/)

## Requirements

- Perl with thread support (check with `perl -V:useithreads`, it should print `useithreads='define'`)
- `g++` and `make`
- Python 3 with numpy, scipy and scikit-learn (only for the evaluation script, see [Evaluation of GREMA](#evaluation-of-grema))

## Getting started

```shell
git clone https://github.com/NctuICLab/GREMA.git
cd GREMA
cd EMA_HFODE; make; cd ..
```

To use the S-system model (`-m s-system`), also compile it:

```shell
cd EMA_SSYSTEM; make; cd ..
```

## Input Data Format

GREMA needs 2 tab-separated input files: 1) time-series data and 2) domain knowledge.

### Time-series data

Example: [Dream4 insilico_size10_1](input/Dream4_10_1_timeseries_expression.txt)

```
repeat_number=5
timepoint_number=21
rep1<TAB>G1<TAB>0.6665<TAB>0.3258<TAB>...(one value per time point)
rep1<TAB>G2<TAB>0.1272<TAB>0.1218<TAB>...
...
rep5<TAB>G10<TAB>...
```

- `repeat_number=`: number of replicates (time series).
- `timepoint_number=`: number of time points in each time series.
- Then one line per replicate and gene: `rep<N>`, gene name, and the expression value at each time point. Every gene must appear in every replicate, and each line must have exactly `timepoint_number` values.

### Domain knowledge

Examples: [no prior knowledge](input/insilico_size10_no_prior_knowledge.txt), [partial knowledge](input/insilico_size10_1_know_knowledge.txt)

```
TF<TAB>GENE<TAB>REGULATORY
G1<TAB>G1<TAB>0
G2<TAB>G1<TAB>?
G3<TAB>G1<TAB>?
...
```

There is one line for every TF-gene pair (including `G1 -> G1`). The `REGULATORY` column can be:

| Value | Meaning |
|-------|---------|
| `?` | unknown, let GREMA infer it |
| `+` | known activation |
| `-` | known repression |
| `0` | known to have no regulation |

If you have no prior knowledge, use the "no prior knowledge" example as a template: all pairs are `?` except the self-regulations (`G1 -> G1`, ...) which are `0`. The "partial knowledge" example additionally sets some pairs to `0`, which makes the search easier.

## Usage of GREMA

```shell
perl GREMA_main.pl -h
Usage: GREMA_main.pl [Options]
Options:
       -i      [FILE]  time-series profile
       -o      [PATH]  Output Directory
       -k      [FILE]  knowledge of regulatory
       -m      [model] {s-system,HFODE} Default is HFODE
       -f      [Fitness] {0,1,2,3,4} Default is 1
       -c      [cc]    {0,1} Default is 0, consider correlation coefficient
       -g      [Generation] Number of generation, default is 1000
       -t      [No]    Number of threads
       -min    [value] ex: 0.01 [min value] default is 0
       -h      Show the usage
```

| Option | Required | Description |
|--------|----------|-------------|
| `-i` | yes | Time-series file. |
| `-k` | yes | Domain knowledge file. |
| `-o` | yes | Output directory. **It must already exist** (GREMA does not create it) and must not contain a `final_results.txt` yet. |
| `-t` | yes | Number of threads, i.e. how many genes are optimised at the same time. Use at most the number of CPU cores. |
| `-m` | no | Model, see [Models](#models). Default `HFODE`. |
| `-g` | no | Maximum number of generations of the genetic algorithm for each run, at least 100. Default 1000. A run also stops early when the best fitness has not changed for 1000 generations. |
| `-f` | no | Fitness function (error between the observed and the simulated expression): `0` squared relative error, `1` squared error (default), `2` squared error of the derivative + squared relative error, `3` squared error of the derivative + squared error, `4` squared error of the derivative only. |
| `-c` | no | `1` adds the correlation-coefficient term `1 - CC` to the fitness. Default `0`. |
| `-min` | no | A parameter whose absolute value is not larger than this value is treated as "no regulation" when determining the regulations. Default `0`. |

## Models

| Model | Option | Executable | Notes |
|-------|--------|------------|-------|
| Hill function-based ODE | `-m HFODE` (default) | `EMA_HFODE/EMA_HFODE` | Regulations are determined from the sign of the Hill coefficients `n_ij`. |
| S-system | `-m s-system` | `EMA_SSYSTEM/EMA_ssystem` | Regulations are determined from the sign of the kinetic orders `g_ij`. The `h_ij` are fixed to 0 (the degradation term is a constant `beta`), and at most `N/3` regulators are searched for each gene (`N` is the number of genes). |

## Run GREMA

### Quick test

Create a new output directory and run a short test with a small number of generations (`-g 100`):

```shell
mkdir -p output/test_run
perl GREMA_main.pl \
    -i input/Dream4_10_1_timeseries_expression.txt \
    -k input/insilico_size10_no_prior_knowledge.txt \
    -o output/test_run/ \
    -g 100 \
    -t 4
```

The progress is printed to the terminal. On a laptop, a 10-gene network with `-g 100` took about 1 minute with HFODE and about 3 minutes with the S-system model.

### Full run

Remove `-g 100` to use the default 1000 generations (the setting used for the results in `output/`). This takes much longer, and the running time grows quickly with the number of genes, so use as many threads as you have cores:

```shell
mkdir -p output/Dream4_10_1_new
perl GREMA_main.pl \
    -i input/Dream4_10_1_timeseries_expression.txt \
    -k input/insilico_size10_no_prior_knowledge.txt \
    -o output/Dream4_10_1_new/ \
    -t 10
```

To use the S-system model, add `-m s-system`.

The `output/Dream4_*` directories in this repository already contain the results of GREMA on the DREAM4 networks. Choose a new output directory, otherwise GREMA stops with `This output folder already has the final_results.txt`.

## Results of GREMA

```shell
less -S output/test_run/final_results.txt
```

Example: [results of Dream4 size10 Network1](output/Dream4_10_1/final_results.txt)

| Column | Description |
|--------|-------------|
| `TF` | Regulator. |
| `GENE` | Target gene. |
| `REGULATORY` | `+` activation, `-` repression, `0` no regulation. |
| `CONFIDENCE_LEVEL` | `1 / iteration`, where *iteration* is the round of GREMA in which the regulation was determined. `1.000` means it was determined in the first round (the most reliable), `0.500` in the second round, `0.333` in the third, and so on. |

The other files in the output directory:

- `G*.txt`: the parameters, fitness and correlation coefficient of the 30 runs of the evolutionary algorithm for each gene (from the final quantification step).
- `tmp/`: intermediate files, such as the knowledge of each round (`*_ForStep*.txt`) and the simulated profiles (`calprofile*`).

## Evaluation of GREMA

Install [scikit-learn](https://scikit-learn.org/stable/install.html) for the evaluation script

### Mac OSX

```shell
pip3 install -U numpy scipy scikit-learn
```

### Linux

```shell
sudo apt-get install build-essential python3-dev python3-setuptools \
                     python3-numpy python3-scipy \
                     libatlas-dev libatlas3gf-base
sudo pip3 install -U scikit-learn
```

To check the usage of evaluation

```shell
python3 evaluation/dream4_evaluate.py --help
usage: Evaluation of GREMA performance. [-h] [-d DIRECTED] [-u UNDIRECTED] input
positional arguments:
  input                 Prediction file (final_results.txt) with path.

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTED, --directed DIRECTED
                        The gold-standard directed GRN file with path
  -u UNDIRECTED, --undirected UNDIRECTED
                        The gold-standard undirected GRN file with path
```

GREMA can reconstruct the **directed GRNs**, so we can evaluate the results using **directed** or **undirected** gold standard GRNs. Here are two types of evaluating commands:

Run the evaluation script for **undirected** gold standard GRN.

```shell
python3 evaluation/dream4_evaluate.py \
    -u evaluation/gold_standards_undirected/10/DREAM4_GoldStandardUndirected_InSilico_Size10_1.tsv \
    output/Dream4_10_1/final_results.txt
==============================
Start evaluating the results of GREMA
Undirected golden file is: evaluation/gold_standards_undirected/10/DREAM4_GoldStandardUndirected_InSilico_Size10_1.tsv
Prediction file is output/Dream4_10_1/final_results.txt
loading gold file:evaluation/gold_standards_undirected/10/DREAM4_GoldStandardUndirected_InSilico_Size10_1.tsv
loading predict results:output/Dream4_10_1/final_results.txt
==============================
Performance of GREMA:
TP=11, TN=23, FP=9, FN=2
Sensitivity=0.846
Specificity=0.719
Accuracy=0.756
ROC AUC=0.782
PR AUC=0.720
==============================

```

Run the evaluation script for **directed** gold standard GRN.

```shell
python3 evaluation/dream4_evaluate.py \
  -d evaluation/gold_standards/10/DREAM4_GoldStandard_InSilico_Size10_1.tsv \
  output/Dream4_10_1/final_results.txt
==============================
Start evaluating the results of GREMA
Directed golden file is: evaluation/gold_standards/10/DREAM4_GoldStandard_InSilico_Size10_1.tsv
Prediction file is output/Dream4_10_1/final_results.txt
loading gold file:evaluation/gold_standards/10/DREAM4_GoldStandard_InSilico_Size10_1.tsv
loading predict results:output/Dream4_10_1/final_results.txt
==============================
Performance of GREMA:
TP=11, TN=57, FP=18, FN=4
Sensitivity=0.733
Specificity=0.760
Accuracy=0.756
ROC AUC=0.747
PR AUC=0.579
==============================
```

## Contact

- Shinn-Ying Ho: syho@nctu.edu.tw
- Ming-Ju Tsai: milutsai.bi98g@g2.nctu.edu.tw
- Jyun-Rong Wang: bengdex.bi99g@g2.nctu.edu.tw
