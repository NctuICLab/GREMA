#!/usr/bin/python3
import argparse
import os
import re
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc
import numpy as np

def natural_key(name):
    """Sort key of a gene name that handles any prefix and any number of digits (G2 < G10)."""
    return tuple((0, int(part), '') if part.isdigit() else (1, 0, part)
                 for part in re.split(r'(\d+)', name) if part != '')


def parse_arguments():
    parser = argparse.ArgumentParser(prog='Evaluation of GREMA performance.',
                                        description='')

    # add long and short argument
    parser.add_argument('input', type=str, nargs=1,
                        help="Prediction file (final_results.txt) with path.")
    parser.add_argument('-d', '--directed', type=str, help='The gold-standard directed GRN file with path')
    parser.add_argument('-u', '--undirected', type=str, help='The gold-standard undirected GRN file with path')
    parser.add_argument('--score', choices=['binary', 'confidence', 'vote'], default='binary',
                        help='Score used for ROC AUC / PR AUC: binary (default, the predicted 0/1 regulation), '
                             'confidence (CONFIDENCE_LEVEL if regulated else 0) or vote (VOTE_SCORE, 5th column). '
                             'TP/TN/FP/FN, sensitivity, specificity and accuracy are always binary.')
    # read arguments from the command line
    args = parser.parse_args()
    if not args.directed and not args.undirected:
        parser.error('one of -d/--directed or -u/--undirected (gold-standard file) is required')
    if args.directed and args.undirected:
        parser.error('use only one of -d/--directed and -u/--undirected')
    return args


def ratio(numerator, denominator):
    """Return numerator/denominator, or None when it is undefined (empty class)."""
    return numerator / denominator if denominator else None


def fmt(value):
    return 'NA' if value is None else '%.3f' % value


def main():
    print('='*30)
    print('Start evaluating the results of GREMA')
    args = parse_arguments()
    prediction_value = {}   # key -> binary prediction
    prediction_score = {}   # key -> score used for the AUCs
    ans_value = {}
    if args.directed:
        print(f'Directed golden file is: {args.directed}')
        evaluate_type = 'directed'
        gold = args.directed
    else:
        print(f'Undirected golden file is: {args.undirected}')
        evaluate_type = 'undirected'
        gold = args.undirected
    if os.path.exists(args.input[0]):
        predict = args.input[0]
        print(f'Prediction file is {predict}')
    else:
        raise SystemExit('input file does not exist')
    if not os.path.exists(gold):
        raise SystemExit(f'gold file does not exist: {gold}')

    def make_key(tf, gene):
        # undirected: the pair is stored in the natural order of the gene names
        if evaluate_type == 'undirected' and natural_key(tf) > natural_key(gene):
            tf, gene = gene, tf
        return (tf, gene)

    print(f'loading gold file:{gold}')
    with open(gold, 'r') as rgold:
        for line_no, line in enumerate(rgold, 1):
            clean_line = line.rstrip('\r\n')
            if not clean_line.strip():
                continue
            fields = clean_line.split('\t')
            if len(fields) < 3:
                raise SystemExit(f'gold file {gold} line {line_no}: expected TF<tab>GENE<tab>0/1')
            tf = fields[0]
            gene = fields[1]
            try:
                regulation = int(fields[2])
            except ValueError:
                raise SystemExit(f'gold file {gold} line {line_no}: regulation must be 0 or 1, got {fields[2]!r}')
            if regulation not in (0, 1):
                raise SystemExit(f'gold file {gold} line {line_no}: regulation must be 0 or 1, got {regulation}')
            if tf == gene:
                continue
            key = make_key(tf, gene)
            if key not in ans_value:
                ans_value[key] = regulation
            else:
                ans_value[key] = max(ans_value[key], regulation)

    print(f'loading predict results:{predict}')
    with open(predict, 'r') as rpredict:
        for line_no, line in enumerate(rpredict, 1):
            clean_line = line.rstrip('\r\n')
            if not clean_line.strip():
                continue
            fields = clean_line.split('\t')
            if len(fields) < 3:
                raise SystemExit(f'prediction file {predict} line {line_no}: expected at least 3 tab separated columns')
            if fields[2] == 'REGULATORY':
                continue
            tf = fields[0]
            gene = fields[1]
            if tf == gene:
                continue
            regulation = 0 if fields[2] == '0' else 1
            if args.score == 'binary':
                score = float(regulation)
            elif args.score == 'confidence':
                if len(fields) < 4:
                    raise SystemExit(f'prediction file {predict} line {line_no}: CONFIDENCE_LEVEL (4th column) is missing')
                try:
                    score = float(fields[3]) if regulation else 0.0
                except ValueError:
                    raise SystemExit(f'prediction file {predict} line {line_no}: bad CONFIDENCE_LEVEL {fields[3]!r}')
            else:
                if len(fields) < 5:
                    raise SystemExit(f'prediction file {predict} line {line_no}: VOTE_SCORE (5th column) is missing, '
                                     'was the file produced by a GREMA version that writes it?')
                try:
                    score = float(fields[4])
                except ValueError:
                    raise SystemExit(f'prediction file {predict} line {line_no}: bad VOTE_SCORE {fields[4]!r}')
            if score != score:
                raise SystemExit(f'prediction file {predict} line {line_no}: score is NaN')
            key = make_key(tf, gene)
            if key not in prediction_value:
                prediction_value[key] = regulation
                prediction_score[key] = score
            else:
                # undirected: both directions of a pair are merged, a regulation in either direction counts
                prediction_value[key] = max(prediction_value[key], regulation)
                prediction_score[key] = max(prediction_score[key], score)

    missing = [k for k in prediction_value if k not in ans_value]
    if missing:
        shown = ', '.join(f'{a}->{b}' for a, b in sorted(missing, key=lambda k: (natural_key(k[0]), natural_key(k[1])))[:5])
        raise SystemExit(f'{len(missing)} predicted pair(s) are not in the gold standard {gold} (e.g. {shown}); '
                         'check that the gene names and the directed/undirected choice match')
    not_predicted = len(ans_value) - len(prediction_value)
    if not_predicted:
        print(f'Warning: {not_predicted} gold pair(s) have no prediction and are ignored')
    if not prediction_value:
        raise SystemExit('no predictions to evaluate')

    predictionofTuples = sorted(prediction_value, key=lambda k: (natural_key(k[0]), natural_key(k[1])))
    predict_list = []
    score_list = []
    ans_list = []
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    for key in predictionofTuples:
        if ans_value[key] == 1:
            if prediction_value[key] == 1:
                TP += 1
            else:
                FN += 1
        else:
            if prediction_value[key] == 0:
                TN += 1
            else:
                FP += 1
        predict_list.append(prediction_value[key])
        score_list.append(prediction_score[key])
        ans_list.append(ans_value[key])

    sen = ratio(TP, TP + FN)
    spe = ratio(TN, TN + FP)
    acc = ratio(TN + TP, TN + TP + FN + FP)

    ans_np = np.array(ans_list)
    score_np = np.array(score_list)
    if len(set(ans_list)) < 2:
        # ROC/PR curves need both positives and negatives in the gold standard
        roc_auc = None
        pr_auc = None
    else:
        roc_auc = roc_auc_score(ans_np, score_np)
        precision, recall, _ = precision_recall_curve(ans_np, score_np)
        pr_auc = auc(recall, precision)
    print('='*30)
    print('Performance of GREMA:')
    print(f'TP={TP}, TN={TN}, FP={FP}, FN={FN}')
    print('Sensitivity=%s' % fmt(sen))
    print('Specificity=%s' % fmt(spe))
    print('Accuracy=%s' % fmt(acc))
    print('ROC AUC=%s' % fmt(roc_auc))
    print('PR AUC=%s' % fmt(pr_auc))
    if args.score != 'binary':
        print(f'(AUCs computed from the {args.score} score)')
    if roc_auc is None:
        print('Note: the gold standard has only one class, AUCs are undefined (NA)')
    print('='*30)

if __name__ == "__main__":
    main()
