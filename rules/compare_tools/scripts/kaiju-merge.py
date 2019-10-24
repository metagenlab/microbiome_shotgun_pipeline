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
    tab=tab.replace(np.nan,0)
    taxid_list = tab['taxid']
    percents = tab['read_percent']
    counts = tab['reads_assigned']
    count_dic = dict(zip(taxid_list, counts))
    percent_dic = dict(zip(taxid_list, percents))
    for taxid in taxid_list:
        tax[taxid] = {'taxid':str(taxid), 'read_percent': percent_dic[taxid], 'read_counts': count_dic[taxid]}
        if taxid > 0:
            scientific_name = ncbi.translate_to_names([taxid])[0]
            tax[taxid]['scientific_name']=scientific_name
            lineage = ncbi.get_lineage(taxid)
            names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)
            path = []
            for sub_taxid in lineage:
                rank = ranks[sub_taxid]
                name = names[sub_taxid]
                if rank in target_ranks:
                    tax[taxid][rank] = name
                    path.append(str(sub_taxid))
            last_taxid = path[len(path) - 1]
            last_rank = ranks[int(last_taxid)]
            tax[taxid]['last_taxid'] = last_taxid
            tax[taxid]['last_rank'] = last_rank
            tax[taxid]['taxid_path'] = '|'.join(path)
        else:
            for j in target_ranks:
                tax[taxid][j] = np.nan
                tax[taxid]['taxid_path'] = np.nan
                tax[taxid]['scientific_name']=np.nan
                tax[taxid]['last_taxid'] = np.nan
                tax[taxid]['last_rank'] = np.nan
                tax[taxid]['taxid_path'] = np.nan
    df = pd.DataFrame.from_dict(tax, orient='index')
    return df


target_ranks = ['superkingdom','phylum','order','family','genus','species']


tables_list=[]
for sample in dic.keys():
    sample_tb=pd.read_csv(dic[sample],sep='\t')
    sample_tb.columns=["file", "read_percent", "reads_assigned", "taxid", "lineage"]
    lin_tax_tb=get_lin_tax(sample_tb,ncbi,target_ranks)
    lin_tax_tb=lin_tax_tb.rename(columns={'read_counts': f'{sample}_counts','read_percent':f'{sample}_percent'})
    lin_tax_tb=lin_tax_tb.groupby(['taxid','scientific_name', 'species', 'genus', 'family', 'order', 'phylum', 'superkingdom', 'taxid_path','last_taxid','last_rank']).sum()
    lin_tax_tb.to_csv(f"benchmark_tools/tables/kaiju/{sample}.tsv",sep='\t')
    tables_list.append(lin_tax_tb)

concatDf=pd.concat(tables_list,sort=False,axis=1)
concatDf.to_csv(snakemake.output[0],sep='\t')