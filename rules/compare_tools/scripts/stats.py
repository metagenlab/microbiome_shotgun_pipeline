import numpy as np
import pandas as pd

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',filename=snakemake.log[0], level=logging.DEBUG)

rank=snakemake.config["target_rank"]
tool_name=snakemake.wildcards.tool


def get_precision_recall_f1(true_tb, tool_tb, rank, tool_name, superkingdom):
    tool_tb = tool_tb.replace(np.nan, 'NA')
    true_tb = true_tb.replace(np.nan, 'NA')

    tool_tb = tool_tb[tool_tb.superkingdom == superkingdom]
    true_tb = true_tb[true_tb.superkingdom == superkingdom]
    names = list(true_tb[f'{rank}'])
    taxids = list(true_tb[f'{rank}_taxid'])
    taxid2names = dict(zip(taxids, names))
    true_list = list(true_tb.groupby([f'{rank}_taxid'], as_index=False).sum()[f'{rank}_taxid'])
    tool_list = list(tool_tb.groupby([f'{rank}_taxid'], as_index=False).sum()[f'{rank}_taxid'])

    tp_list = [i for i in true_list if i in true_list and i in tool_list]
    fp_list = [i for i in tool_list if i not in true_list and i in tool_list]
    fn_list = [i for i in true_list if i in true_list and i not in tool_list]
    tp = len(tp_list)
    fp = len(fp_list)
    fn = len(fn_list)
    try:
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1 = 2 * (precision * recall) / (precision + recall)
    except ZeroDivisionError:
        precision = 0
        recall = 0
        f1 = 0
    scores = {}
    scores[superkingdom] = {'tool': tool_name, 'tp': tp, 'fp': fp, 'fn': fn, 'precision': precision, 'recall': recall,
                            'F1_score': f1}
    if fn > 0:
        for false_neg in fn_list:
            logging.info(f'false negative found for {superkingdom}: {taxid2names[false_neg]}')
    df = pd.DataFrame.from_dict(scores, orient='index')
    return df


def gold_standard_confusion_table(true_tb,tool_tb,rank,tool_name,superkingdom):
    tool_tb = tool_tb.replace(np.nan, 'NA')
    true_tb = true_tb.replace(np.nan, 'NA')

    tool_tb = tool_tb[tool_tb.superkingdom == superkingdom]
    true_tb = true_tb[true_tb.superkingdom == superkingdom]
    names = list(true_tb[f'{rank}'])
    taxids = list(true_tb[f'{rank}_taxid'])
    taxid2names = dict(zip(taxids, names))
    true_list = list(true_tb.groupby([f'{rank}_taxid'], as_index=False).sum()[f'{rank}_taxid'])
    tool_list = list(tool_tb.groupby([f'{rank}_taxid'], as_index=False).sum()[f'{rank}_taxid'])
    tp_list = [i for i in true_list if i in true_list and i in tool_list]
    fn_list = [i for i in true_list if i in true_list and i not in tool_list]

    return

gold_standard=pd.read_csv(snakemake.input.gold_standard,sep='\t')
tool_output=pd.read_csv(snakemake.input.tool_out,sep='\t')

virus_scores=get_precision_recall_f1(gold_standard,tool_output,rank,tool_name,'Viruses')
bacteria_scores=get_precision_recall_f1(gold_standard,tool_output,rank,tool_name,'Bacteria')
scores_table=pd.concat([virus_scores,bacteria_scores],sort=False,axis=0)
scores_table.to_csv(snakemake.output[0],sep='\t')