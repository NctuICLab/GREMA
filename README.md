# GREMA

**GREMA** (**G**ene networks **R**econstruction using **E**volutionary **M**odelling **A**lgorithm) is a program for inferring  a novel type of gene regulatory network (GRN) with confidence levels for every inferred regulation, which is emulated GRN (eGRN). The higher the confidence level, the more accurate the inferred regulation. GREMA gradually determines the regulations of an eGRN with confidence levels in descending order using either an S-system or a Hill function-based ordinary differential equation model. It makes use of an evolutionary modelling algorithm (EMA) that is based on evolutionary intelligence, including both crowd wisdom and an evolutionary strategy, to cope with the underdetermined problem. EMA uses an intelligent genetic algorithm to solve the large-scale parameter optimisation problem.

## Citing GREMA

GREMA: Modelling of emulated gene regulatory networks with confidence levels based on evolutionary intelligence to cope with the underdetermined problem. _Bioinformatics_, 2020, [https://doi.org/10.1093/bioinformatics/btaa267](https://doi.org/10.1093/bioinformatics/btaa267)

## Input Data Format

There are 2 input tsv files: 1) Time-series data and 2) Domain knowledge data.

**The format of time-series data** (e.g., [Dream4 insilico_size10_1](input/Dream4_10_1_timeseries_expression.txt))

- repeat_number= (e.g.,repeat_number=5)
- timepoint_number= (e.g.,timepoint_number=21)
- time-series profile (repeat_number \<tab\> gene_name \<tab\> gene_expression)

**The format of domain knowledge data** (e.g., [Dream 4 insilico size10_1](input/insilico_size10_1_know_knowledge.txt))

- TF
- GENE
- REGULATORY (+:activation, -:repression, ?:unknown, 0:No regulation).

## Getting start

 ```shell
 git clone https://github.com/NctuICLab/GREMA.git
 cd GREMA
 cd EMA_HFODE;make
 ```

## Usage of GREMA

To check the usage of GREMA_main.pl

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

## Run GREMA

```shell
perl GREMA_main.pl \
    -i input/Dream4_10_1_timeseries_expression.txt \
    -k input/insilico_size10_no_prior_knowledge.txt \
    -o output/Dream4_10_1/ \
    -t 10
```

## Results of GREMA

```shell
less -S output/Dream4_10_1/final_results.txt
```

The format of results of GREMA (e.g., Results of Dream4 size10 Network1: [final results](output/Dream4_10_1/final_results.txt))

- TF
- GENE
- REGULATORY (+:activation, -:repression, ?:unknown, 0:No regulation).
- CONFIDENCE_LEVEL

## Evaluation of GREMA

Installing [Scikit learn](https://scikit-learn.org/0.16/install.html) for evaluation script

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
    -u evalutation/gold_standards_undirected/10/DREAM4_GoldStandardUndirected_InSilico_Size10_1.tsv \
    output/Dream4_10_1/final_results.txt
==============================
Start evaluating the results of GREMA
Undirected golden file is: evalutation/gold_standards_undirected/10/DREAM4_GoldStandardUndirected_InSilico_Size10_1.tsv
Prediction file is output/Dream4_10_1/final_results.txt
loading gold file:evalutation/gold_standards_undirected/10/DREAM4_GoldStandardUndirected_InSilico_Size10_1.tsv
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
  -d evalutation/gold_standards/10/DREAM4_GoldStandard_InSilico_Size10_1.tsv \
  output/Dream4_10_1/final_results.txt
==============================
Start evaluating the results of GREMA
Directed golden file is: evalutation/gold_standards/10/DREAM4_GoldStandard_InSilico_Size10_1.tsv
Prediction file is output/Dream4_10_1/final_results.txt
loading gold file:evalutation/gold_standards/10/DREAM4_GoldStandard_InSilico_Size10_1.tsv
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
