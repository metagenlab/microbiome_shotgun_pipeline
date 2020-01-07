import numpy as np
import pandas as pd

rank=snakemake.config["target_rank"]
tool_name=snakemake.wildcards.tool

threshold=snakemake.params.threshold
superkingdom=snakemake.params.superkingdom

def get_precision_recall_f1(true_tb, tool_tb, rank, tool_name):
    if rank=='scientific_name':
        taxid_col='taxid'
    else:
        taxid_col=f'{rank}_taxid'
    true_taxids = list(set(true_tb[taxid_col]))
    tool_taxids= list(set(tool_tb[taxid_col]))

    true_list=[str(int(float(taxid))) for taxid in true_taxids]
    tool_list = [str(int(float(taxid))) for taxid in tool_taxids]

    tp=0
    fn=0
    for i in true_list:
        if tp < len(true_list) and i in true_list and i in tool_list:
            tp += 1
        if tp < len(true_list) and i in true_list and i not in tool_list:
            fn += 1
        if tp + fn == len(true_list):
            break
    fp = len(tool_list) - tp
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

    true_taxids = list(true_tb[taxid_col])
    true_names = list(true_tb[f'{rank}'])
    tool_taxids=list(tool_tb[taxid_col])
    tool_names=list(tool_tb[f'{rank}'])
    tool_samples=list(tool_tb['sample'])
    true_samples=list(true_tb['sample'])
    true_list = [str(int(float(taxid))) for taxid in true_taxids]
    tool_list = [str(int(float(taxid))) for taxid in tool_taxids]#some tools have taxids as floats, this transforms them to strings of integers
    true_taxid_names = dict(zip(true_list, true_names))
    tool_taxid_names = dict(zip(tool_list, tool_names))
    tool_taxid_samples=dict(zip(tool_list,tool_samples))
    true_taxid_samples=dict(zip(true_list,true_samples))
    tp_list = [i for i in true_list if i in true_list and i in tool_list]
    fp_list = [i for i in tool_list if i not in true_list and i in tool_list]
    fn_list = [i for i in true_list if i in true_list and i not in tool_list]

    dic={}
    for taxid in true_taxid_names:
        dic[taxid]={}
        if taxid in tp_list:
            dic[taxid]={'tool':tool_name,'tool-sample':tool_taxid_samples[taxid],'true-sample':true_taxid_samples[taxid],'taxid':str(taxid),f'{rank}':true_taxid_names[taxid],'present':"TP"}
        if taxid in fn_list:
            dic[taxid]={'tool':tool_name,'tool-sample':'NA','true-sample':true_taxid_samples[taxid],'taxid':str(taxid),f'{rank}':true_taxid_names[taxid],'present':"FN"}
    for taxid in fp_list:
            dic[taxid]={'tool':tool_name,'tool-sample':tool_taxid_samples[taxid],'true-sample':'NA','taxid':str(taxid),f'{rank}':tool_taxid_names[taxid],'present':"FP"}

    df=pd.DataFrame.from_dict(dic, orient='index')
    return df

dtype={'taxid':'object','superkingdom_taxid':'object','phylum_taxid':'object','order_taxid':'object',
       'family_taxid':'object','genus_taxid':'object','species_taxid':'object'}#Need to set all taxids to object, otherwise groupby will sum them
gold_standard=pd.read_csv(snakemake.input.gold_standard,sep='\t',dtype=dtype)
tool_output=pd.read_csv(snakemake.input.tool_out,sep='\t',dtype=dtype)

true_superkingdom=gold_standard[gold_standard.superkingdom==superkingdom]
true_superkingdom=true_superkingdom.groupby(['sample',f'{rank}',f'{rank}_taxid'],as_index=False).sum()
true_superkingdom = true_superkingdom.replace(np.nan, 'NA')

tool_superkingdom=tool_output[tool_output.superkingdom==superkingdom]

tool_superkingdom.insert(loc=tool_superkingdom.shape[1],column='tool',value=[tool_name]*len(tool_superkingdom))
tool_superkingdom=tool_superkingdom.groupby(['tool','sample',f'{rank}_taxid',f'{rank}'],as_index=False).sum()
tool_superkingdom=tool_superkingdom[tool_superkingdom.read_counts>=threshold]
tool_superkingdom=tool_superkingdom.replace(np.nan,'NA')



scores=get_precision_recall_f1(true_superkingdom,tool_superkingdom,rank,tool_name)
presence=get_presence(true_superkingdom,tool_superkingdom,rank,tool_name)
scores.to_csv(snakemake.output.scores,sep='\t',index=None)
presence.to_csv(snakemake.output.presence,sep='\t',index=None)
