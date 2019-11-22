import numpy as np
import pandas as pd

rank=snakemake.config["target_rank"]
tool_name=snakemake.wildcards.tool

bac_threshold=snakemake.params.bac_threshold
vir_threshold=snakemake.params.vir_threshold

def get_precision_recall_f1(true_tb, tool_tb, rank, tool_name):
    if rank=='scientific_name':
        taxid_col='taxid'
    else:
        taxid_col=f'{rank}_taxid'
    true_list = list(true_tb.groupby([taxid_col], as_index=False).sum()[taxid_col])
    tool_list = list(tool_tb.groupby([taxid_col], as_index=False).sum()[taxid_col])
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
    scores[tool_name] = {'tool':tool_name,'tp': tp, 'fp': fp, 'fn': fn, 'precision': precision, 'recall': recall,
                            'F1_score': f1}
    df = pd.DataFrame.from_dict(scores, orient='index')
    return df



def get_presence(true_tb,tool_tb,rank,tool_name):

    if rank=='scientific_name':
        taxid_col='taxid'
    else:
        taxid_col=f'{rank}_taxid'

    true_names = list(true_tb[f'{rank}'])
    true_taxids = list(true_tb[taxid_col])
    true_taxid_names = dict(zip(true_taxids, true_names))
    tool_taxids=list(tool_tb[taxid_col])
    tool_names=list(tool_tb[f'{rank}'])
    tool_taxid_names=dict(zip(tool_taxids,tool_names))
    true_list = list(true_tb.groupby([taxid_col], as_index=False).sum()[taxid_col])
    tool_list = list(tool_tb.groupby([taxid_col], as_index=False).sum()[taxid_col])
    tp_list = [i for i in true_list if i in true_list and i in tool_list]
    fn_list = [i for i in true_list if i in true_list and i not in tool_list]
    fp_list = [i for i in tool_list if i not in true_list and i in tool_list]
    dic={}
    for taxid in true_taxid_names:
        dic[taxid]={}
        if taxid in tp_list:
            dic[taxid]={'tool':tool_name,'taxid':str(taxid),f'{rank}':true_taxid_names[taxid],'present':"TP"}
        if taxid in fn_list:
            dic[taxid]={'tool':tool_name,'taxid':str(taxid),f'{rank}':true_taxid_names[taxid],'present':"FN"}
    for taxid in fp_list:
            dic[taxid]={'tool':tool_name,'taxid':str(taxid),f'{rank}':tool_taxid_names[taxid],'present':"FP"}

    df=pd.DataFrame.from_dict(dic, orient='index')
    return df

gold_standard=pd.read_csv(snakemake.input.gold_standard,sep='\t')
tool_output=pd.read_csv(snakemake.input.tool_out,sep='\t')

true_vir=gold_standard[gold_standard.superkingdom=='Viruses']
true_bac=gold_standard[gold_standard.superkingdom=='Bacteria']

tool_vir=tool_output[tool_output.superkingdom=='Viruses']
tool_vir=tool_vir[tool_vir.read_counts>vir_threshold]
tool_vir=tool_vir.replace(np.nan,'NA')
tool_bac=tool_output[tool_output.superkingdom=='Bacteria']
tool_bac=tool_bac[tool_bac.read_counts>bac_threshold]
tool_bac=tool_bac.replace(np.nan,'NA')

gold_standard = gold_standard.replace(np.nan, 'NA')

virus_scores=get_precision_recall_f1(true_vir,tool_vir,rank,tool_name)
virus_presence=get_presence(true_vir,tool_vir,rank,tool_name)
virus_scores.to_csv(snakemake.output.vir_scores,sep='\t',index=None)
virus_presence.to_csv(snakemake.output.vir_presence,sep='\t',index=None)
bacteria_scores=get_precision_recall_f1(true_bac,tool_bac,rank,tool_name)
bacteria_presence=get_presence(true_bac,tool_bac,rank,tool_name)
bacteria_presence.to_csv(snakemake.output.bac_presence,sep='\t',index=None)
bacteria_scores.to_csv(snakemake.output.bac_scores,sep='\t',index=None)