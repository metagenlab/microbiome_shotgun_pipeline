import pandas as pd
from ete3 import NCBITaxa
import numpy as np
ncbi = NCBITaxa()


input_file_list=snakemake.input
dic={}
for path in input_file_list:
    sample=path.split('/')[1]
    dic[sample]=path


def get_lin_tax(tab, ncbi, target_ranks):
    tax = {}
    taxid_list = tab['taxid']
    percents = tab['read_percent']
    counts = tab['reads_assigned']
    count_dic = dict(zip(taxid_list, counts))
    percent_dic = dict(zip(taxid_list, percents))
    for taxid in taxid_list:
        tax[taxid] = {'taxid': str(taxid), 'read_percent': percent_dic[taxid], 'read_counts': count_dic[taxid]}
        if taxid > 0:
            lin_txid = ncbi.get_lineage(taxid)
            lin_translation = ncbi.get_taxid_translator(lin_txid)
            taxid2rank = ncbi.get_rank(list(lin_translation.keys()))
            taxid_path = []
            for n, k in enumerate(taxid2rank.keys()):
                if taxid2rank[k] in target_ranks:
                    tax[taxid][taxid2rank[k]] = lin_translation[k]
                    tax[taxid]['taxid_path'] = taxid_path.append(str(lin_txid[n]))
            tax[taxid]['taxid_path'] = '|'.join(taxid_path)
        else:
            for j in target_ranks:
                tax[taxid][j] = np.nan
                tax[taxid]['taxid_path'] = np.nan

    df = pd.DataFrame.from_dict(tax, orient='index')
    return df


target_ranks = ['superkingdom','phylum','order','family','genus','species']

tables_list = []

for sample in dic.keys():
    sample_tb=pd.read_csv(dic[sample],sep='\t',names=["read_percent","rooted_frag_num","reads_assigned","rank_code","taxid","scientific_names"])
    lin_tax_tb=get_lin_tax(sample_tb,ncbi,target_ranks)
    lin_tax_tb.rename(columns={'read_counts': f'{sample}_counts','read_percent': f'{sample}_percent'})
    lin_tax_tb.groupby(['taxid', 'species', 'genus', 'family', 'order', 'phylum', 'superkingdom', 'taxid_path'],
                       as_index=False).sum()
    tables_list.append(lin_tax_tb)

concatDf=pd.concat(tables_list,sort=False,axis=0,)
concatDf.to_csv(snakemake.output[0],sep='\t')