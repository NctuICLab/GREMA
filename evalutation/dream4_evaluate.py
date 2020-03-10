#!/usr/bin/python3
import argparse
import os
import re
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(prog='Evaluation of GREMA performance.',
                                        description='')

    # add long and short argument
    parser.add_argument('input', type=str, nargs=1,
                        help="Prediction file (final_results.txt) with path.")
    parser.add_argument('-d', '--directed', type=str, help='The gold-standard directed GRN file with path')
    parser.add_argument('-u', '--undirected', type=str, help='The gold-standard undirected GRN file with path')
    # read arguments from the command line
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    prediction = {}
    prediction_value = {}
    ans = {}
    ans_value = {}
    if args.directed:
        print(f'directed golden file is: {args.directed}')
        evaluate_type = 'directed'
    elif args.undirected:
        print(f'undirected golden file is: {args.undirected}')
        evaluate_type = 'undirected'
    if os.path.exists(args.input[0]):
        predict = args.input[0]
        print(f'input file is {predict}')
    else:
        raise Exception('input file does not exist')

    if evaluate_type == 'undirected':
        gold = args.undirected
    else:
        gold = args.directed
    print(f'load gold file:{gold}')
    rgold = open(gold,'r')
    for line in rgold:
        clean_line = line.rstrip('\r\n')
        fields = clean_line.split('\t')
        tf = fields[0]
        gene = fields[1]
        regulation = int(fields[2])
        if tf == gene:
            continue
        else:
            key = key = tf + "\t" + gene
        if key not in ans:
            ans_value[key] = regulation

    print(f'load predict results:{predict}')
    rpredict = open(predict,'r')
    count = 0
    for line in rpredict:
        clean_line = line.rstrip('\r\n')
        fields = clean_line.split('\t')
        if fields[2] == 'REGULATORY':
            continue
        tf = fields[0]
        tf_number = int(tf.replace('G', ''))
        #print(f"tf_number:{tf_number},tf:{tf}")
        gene = fields[1]
        gene_number = int(gene.replace('G', ''))
        if fields[2] == '0':
            regulation = 0
        else:
            regulation = 1
        #print(f"gene_number:{gene_number},tf:{gene}")
        if evaluate_type == 'undirected':
            if tf_number == gene_number:
                continue
            elif tf_number > gene_number:
                key = gene + "\t" + tf
                value = gene_number * 10 + tf_number
                #origin_key = tf + "\t" + gene
                #print(f"origin key is {origin_key}, new key is: {key}")
            else:
                key = tf + "\t" + gene
                value = tf_number * 10 + gene_number
                #print(f'key is {key}')
        else:
            if tf_number == gene_number:
                continue
            else:
                key = tf + "\t" + gene
                value = tf_number * 10 + gene_number
        if key not in prediction or prediction_value[key] == 0:
            prediction[key] = value
            prediction_value[key] = regulation
            #print(f'prediction_value[{key}]={regulation}')
            count += 1
    
    #print(f'total count:{count}')
    predictionofTuples = sorted(prediction.items(), key=lambda x: x[1])
    i = 0
    predict_list = []
    ans_list = []
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    for elem in predictionofTuples:
        #print(f'{elem[0]} predict:{prediction_value[elem[0]]}, ans:{ans_value[elem[0]]}')
        if ans_value[elem[0]] == 1:
            if prediction_value[elem[0]] == 1:
                TP += 1
            else:
                FN += 1
        elif ans_value[elem[0]] == 0:
            if prediction_value[elem[0]] == 0:
                TN += 1
            else:
                FP += 1
        predict_list.append(prediction_value[elem[0]])
        ans_list.append(ans_value[elem[0]])
        i += 1
    
    sen = TP/(TP + FN)
    spe = TN/(TN + FP)
    acc = (TN + TP)/(TN + TP + FN + FP)
    print(f'TP={TP}, TN={TN}, FP={FP}, FN={FN}')
    print('Sensitivity=%.3f' % (sen))
    print('Specificity=%.3f' % (spe))
    print('Accuracy=%.3f' % (acc))
    
    #print(predict_list)
    #print(ans_list)
    predict_np = np.array(predict_list)
    ans_np = np.array(ans_list)
    #print(predict_np)
    #print(ans_np)
    roc_auc = roc_auc_score(ans_np, predict_np)
    precision, recall, _ = precision_recall_curve(ans_np, predict_np)
    pr_auc = auc(recall, precision)
    print('ROC AUC=%.3f' % (roc_auc))
    print('PR AUC=%.3f' % (pr_auc))

if __name__ == "__main__":
    main()   
