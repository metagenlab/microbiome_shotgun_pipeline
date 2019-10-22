import pandas as pd
from ete3 import NCBITaxa
import numpy as np
ncbi = NCBITaxa()



virus_input=snakemake.input['vir']
bacteria_input=snakemake.input['bac']

vir_dic={}
for path in virus_input:
    sample=path.split('/')[1]
    vir_dic[sample]=path

bac_dic={}
for path in bacteria_input:
    sample=path.split('/')[1]
    bac_dic[sample]=path


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

tables_list=[]
target_ranks = ['superkingdom','phylum','order','family','genus','species']

for sample in vir_dic.keys():
    vir_tb_input=pd.read_csv(vir_dic[sample],sep='\t',names=['taxid','reads_assigned','read_percent','reads_assigned_total','reads_uniquely_assigned','rank','name','NA'])
    vir_tb_linear_tax=get_lin_tax(vir_tb_input,ncbi,target_ranks)
    bac_tb_input=pd.read_csv(bac_dic[sample],sep='\t',names=['taxid','reads_assigned','read_percent','reads_assigned_total','reads_uniquely_assigned','rank','name','NA'])
    bac_tb_linear_tax=get_lin_tax(bac_tb_input,ncbi,target_ranks)
    all_tb=pd.concat([bac_tb_linear_tax, vir_tb_linear_tax], sort=False)
    all_grouped=all_tb.groupby(['taxid','species', 'genus', 'family', 'order', 'phylum', 'superkingdom','taxid_path'],as_index=False).sum()
    all_grouped=all_grouped.rename(columns={'read_counts':f'{sample}_counts','read_percent': f'{sample}_percent'})
    tables_list.append(all_grouped)

all_ganon=pd.concat(tables_list,sort=False,axis=0)

all_ganon.to_csv(snakemake.output[0],sep='\t')