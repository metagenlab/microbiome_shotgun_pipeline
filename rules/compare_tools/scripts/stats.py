import numpy as np
import pandas as pd

rank=snakemake.config["target_rank"]
tool_name=snakemake.wildcards.tool

threshold=snakemake.params.threshold
superkingdom=snakemake.params.superkingdom


def get_presence_per_sample(true_superkingdom, tool_superkingdom, tool_name, sample, rank):

    target_col = f'{rank}_taxid'
    if rank=='scientific_name':
        target_col='taxid'

    if len(list(set(true_superkingdom['sample']))) == 1:  # For the case of my simulated mock community all samples contain the same genomes
        true_taxids = list(true_superkingdom[target_col])
        true_names = list(true_superkingdom[f'{rank}'])
    else:
        true_taxids = list(true_superkingdom[true_superkingdom['sample'] == sample][target_col])
        true_names = list(true_superkingdom[true_superkingdom['sample'] == sample][f'{rank}'])
    true_dic = dict(zip(true_taxids, true_names))

    tool_taxids = list(tool_superkingdom[tool_superkingdom['sample'] == sample][target_col])
    tool_taxids = [str(int(float(taxid))) for taxid in tool_taxids]#for converting float taxid to integers
    tool_names = list(tool_superkingdom[tool_superkingdom['sample'] == sample][f'{rank}'])
    tool_dic = dict(zip(tool_taxids, tool_names))
    dic = {}
    for taxid in list(set(true_taxids)):
        dic[taxid] = {}
        if taxid in tool_taxids:
            dic[taxid] = {'tool': tool_name, 'sample': sample, f'{rank}': true_dic[taxid], 'taxid': taxid,
                          'present': 'TP'}
        if taxid not in tool_taxids:
            dic[taxid] = {'tool': tool_name, 'sample': sample, f'{rank}': true_dic[taxid], 'taxid': taxid,
                          'present': 'FN'}
    for taxid in list(set(tool_taxids)):
        if taxid not in true_taxids:
            dic[taxid] = {'tool': tool_name, 'sample': sample, f'{rank}': tool_dic[taxid], 'taxid': taxid,
                          'present': 'FP'}
    df = pd.DataFrame.from_dict(dic, orient='index')
    return df


def get_precision_recall_f1(tool,tool_tb,threshold):
    matrix={}
    tp=len(tool_tb[tool_tb.present=='TP'])
    fn=len(tool_tb[tool_tb.present=='FN'])
    fp=len(tool_tb[tool_tb.present=='FP'])
    try:
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1 = 2 * (precision * recall) / (precision + recall)
    except ZeroDivisionError:
        precision = 0
        recall = 0
        f1 = 0
    matrix[tool]={'tool':tool,'tp':tp,'fn':fn,'fp':fp,'precision': precision, 'recall': recall,
                            'F1_score': f1,'threshold':threshold}
    mt = pd.DataFrame.from_dict(matrix, orient='index')
    return mt


dtype={'taxid':'object','superkingdom_taxid':'object','phylum_taxid':'object','order_taxid':'object',
       'family_taxid':'object','genus_taxid':'object','species_taxid':'object'}#Need to set all taxids to object, otherwise groupby will sum them
gold_standard=pd.read_csv(snakemake.input.gold_standard,sep='\t',dtype=dtype)
tool_output=pd.read_csv(snakemake.input.tool_out,sep='\t',dtype=dtype)

true_superkingdom=gold_standard[gold_standard.superkingdom==superkingdom]
true_superkingdom=true_superkingdom.groupby(['sample',f'{rank}',f'{rank}_taxid'],as_index=False).sum()
true_superkingdom = true_superkingdom.replace(np.nan, 'NA')



tool_superkingdom=tool_output[tool_output.superkingdom==superkingdom]
tool_superkingdom.insert(loc=tool_superkingdom.shape[1],column='tool',value=[tool_name]*len(tool_superkingdom))


if tool_name=='ezvir':
    tool_superkingdom=tool_superkingdom.groupby(['tool','sample',f'{rank}_taxid',f'{rank}'],as_index=False).sum()
    tool_superkingdom['read_counts']=tool_superkingdom['read_counts']/2
else:
    tool_superkingdom = tool_superkingdom.groupby(['tool', 'sample', f'{rank}_taxid', f'{rank}'], as_index=False).sum()

tool_superkingdom=tool_superkingdom[tool_superkingdom.read_counts>=threshold]
tool_superkingdom=tool_superkingdom.replace(np.nan,'NA')


all_samples=[]
samples=list(set(tool_superkingdom['sample']))
for sample in samples:
    all_samples.append(get_presence_per_sample(true_superkingdom,tool_superkingdom,tool_name,sample,rank))
presence=pd.concat(all_samples)

scores=get_precision_recall_f1(tool_name,presence,threshold)

scores.to_csv(snakemake.output.scores,sep='\t',index=None)
presence.to_csv(snakemake.output.presence,sep='\t',index=None)
