import numpy as np
import pandas as pd

rank='scientific_name'


threshold=snakemake.params.threshold
superkingdom=snakemake.params.superkingdom
tool_name=snakemake.wildcards.tool

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
    true_ct=list(true_superkingdom[true_superkingdom['sample'] == sample]['Ct'])

    true_taxids=[str(int(float(taxid))) for taxid in true_taxids]
    true_dic = dict(zip(true_taxids, true_names))
    taxid2ct = dict(zip(true_taxids, true_ct))

    tool_taxids = list(tool_superkingdom[tool_superkingdom['sample'] == sample][target_col])
    tool_taxids = [str(int(float(taxid))) for taxid in tool_taxids]#for converting float taxid to integers
    tool_names = list(tool_superkingdom[tool_superkingdom['sample'] == sample][f'{rank}'])
    tool_counts = list(tool_superkingdom[tool_superkingdom['sample'] == sample]['read_counts'])

    tool_dic = dict(zip(tool_taxids, tool_names))
    taxid2counts=dict(zip(tool_taxids, tool_counts))
    dic = {}
    for taxid in list(set(true_taxids)):
        dic[taxid] = {}
        if taxid in tool_taxids:
            dic[taxid] = {'tool': tool_name, 'sample': sample, f'{rank}': true_dic[taxid], 'taxid': taxid,
                          'present': 'TP','read_counts':taxid2counts[taxid],'Ct':taxid2ct[taxid]}
        if taxid not in tool_taxids:
            dic[taxid] = {'tool': tool_name, 'sample': sample, f'{rank}': true_dic[taxid], 'taxid': taxid,
                          'present': 'FN','read_counts':'NA','Ct':taxid2ct[taxid]}
    for taxid in list(set(tool_taxids)):
        if taxid not in true_taxids:
            dic[taxid] = {'tool': tool_name, 'sample': sample, f'{rank}': tool_dic[taxid], 'taxid': taxid,
                          'present': 'FP','read_counts':taxid2counts[taxid],'Ct':'NA'}
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
true_superkingdom=gold_standard[gold_standard.superkingdom==superkingdom]
#true_superkingdom=true_superkingdom.groupby(['sample',f'{rank}','taxid'],as_index=False).sum()
true_superkingdom = true_superkingdom.replace(np.nan, 'NA')


tool_out=pd.read_csv(snakemake.input.tool_out,sep='\t')

tool_out['sample']=[sample.split('_')[0] for sample in tool_out['sample']]
tool_out=tool_out[tool_out.superkingdom==superkingdom]
tool_superkingdom = tool_out.groupby(['sample',f'{rank}','taxid'], as_index=False).sum()
tool_out.insert(loc=tool_out.shape[1],column='tool',value=[tool_name]*len(tool_out))

tool_superkingdom.replace(np.nan,'NA',inplace=True)
tool_superkingdom=tool_superkingdom[tool_superkingdom.read_counts>=threshold]




all_samples=[]
samples=list(set(true_superkingdom['sample']))
sample_subset=tool_superkingdom[tool_superkingdom['sample'].isin(samples)]
taxids=list(set(true_superkingdom['taxid']))
sample_taxid_subset=sample_subset[sample_subset['taxid'].isin(taxids)]



for sample in samples:
    all_samples.append(get_presence_per_sample(true_superkingdom,sample_taxid_subset,tool_name,sample,rank))
presence=pd.concat(all_samples)

scores=get_precision_recall_f1(tool_name,presence,threshold)

scores.to_csv(snakemake.output.scores,sep='\t',index=None)
presence.to_csv(snakemake.output.presence,sep='\t',index=None)
