# GREMA
This work proposes an evolutionary modelling algorithm (EMA) that is based on evolutionary intelligence, including both crowd wisdom and an evolutionary strategy, to cope with the underdetermined problem. EMA uses an intelligent genetic algorithm to solve the large-scale parameter optimisation problem. 
An EMA-based method, GREMA, infers a novel type of gene regulatory network with confidence levels for every inferred regulation. The higher the confidence level, the more accurate the inferred regulation. GREMA gradually determines the regulations of an eGRN with confidence levels in descending order using either an S-system or a Hill function-based ordinary differential equation model. 

## Input Data Format
There are 2 input tsv files: time-series data and domain knowledge data.
1. The format of time-series data [Dream4 insilico_size10_1](input/Dream4_10_1_timeseries_expression.txt)
 - repeat_number= (ex:repeat_number=5)
 - timepoint_number= (ex:timepoint_number=21)
 - time-series profile (rep1 \<tab\> gene_name \<tab\> gene_expression)
2. The format of domain knowledge data [Dream 4 insilico size10_1](input/insilico_size10_1_know_knowledge.txt)
 - TF
 - GENE
 - REGULATORY (+:activation, -:repression, ?:unknown, 0:No regulation).
 
 ## Getting start
 ```shell
 git clone https://github.com/NctuICLab/GREMA.git
 cd GREMA
 cd EMA_HFODE
 make
 cd ../
 cd ../
 perl GREMA.main.pl -i input/Dream4_10_1_timeseries_expression.txt -k input/insilico_size10_1_know_knowledge.txt -o output/Dream4_10_1/ -t 10
 ```
 
 ## Usage of GREMA_main.pl
 ```shell
 Usage: GREMA_main.pl [Options]
 Options:
        -i      [FILE]  time-series profile
        -o      [PATH]  Output Directory
        -k      [FILE]  knowledge of regulatory
        -m      [model] {s-system,HFODE}
        -t      [No]    Number of threads
        -h      Show the usage
 ```

## Results of GREMA
```shell
less -S output/Dream4_10_1/final_results.txt
```
The format of results of GREMA [final results](output/final_results.txt)
- TF
- GENE
- REGULATORY (+:activation, -:repression, ?:unknown, 0:No regulation).
- CONFIDENCE_LEVEL

## Evaluation of GREMA
```shell
python3 evalutation/dream4_evaluate.py --help
usage: Evaluation of GREMA performance. [-h] [-d DIRECTED] [-u UNDIRECTED]
                                        input

positional arguments:
  input                 Path for GREMA's prediction file.

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTED, --directed DIRECTED
       Path for gold-standard directed GRN
  -u UNDIRECTED, --undirected UNDIRECTED
       Path for gold-standard undirected GRN
python3 evalutation/dream4_evaluate.py -u evalutation/gold_standards_undirected/10/DREAM4_GoldStandardUndirected_InSilico_Size10_1.tsv output/Dream4_10_1/final_results.txt
oad gold file:evalutation/gold_standards_undirected/10/DREAM4_GoldStandardUndirected_InSilico_Size10_1.tsv
load predict results:output/Dream4_10_1/final_results.txt
TP=11, TN=23, FP=9, FN=2
Sensitivity=0.846
Specificity=0.719
Accuracy=0.756
ROC AUC=0.782
PR AUC=0.720
```




  
