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
            rank_list = list(ranks.values())
            name_list = list(names.values())
            rank2names = dict(zip(rank_list, name_list))

            for sub_taxid in lineage:
                rank = ranks[sub_taxid]
                if rank in target_ranks:
                    tax[taxid][f"{rank}_taxid"] = sub_taxid
            for n, rank in enumerate(target_ranks):
                if rank in rank2names.keys():
                    tax[taxid][rank] = rank2names[rank]
                else:
                    previous_name = 'root'
                    for ranks_to_skip in range(1, len(target_ranks) + 1):
                        previous_rank = target_ranks[n - ranks_to_skip]
                        if previous_rank in rank2names.keys():
                            previous_name = rank2names[previous_rank]
                            break
                        if previous_rank not in rank2names.keys():
                            previous_name = ''
                            continue
                    tax[taxid][rank] = f'{previous_name}_{rank[0:1]}'
    df = pd.DataFrame.from_dict(tax, orient='index')
    df = df.replace(np.nan, 'NA')
    return df

target_ranks = ['superkingdom','phylum','order','family','genus','species']
tables_list=[]
split_path=snakemake.output[0].split('/')
tool_path='/'.join(split_path[0:len(split_path)-1])

for sample in dic.keys():
    sample_tb=pd.read_csv(dic[sample],sep='\t')
    sample_tb.columns=["file", "read_percent", "reads_assigned", "taxid", "lineage"]
    lin_tax_tb=get_lin_tax(sample_tb,ncbi,target_ranks)
    lin_tax_tb['sample']=[sample]*len(lin_tax_tb)
    lin_tax_tb=lin_tax_tb.groupby(['superkingdom','superkingdom_taxid','phylum','phylum_taxid','order','order_taxid','family','family_taxid','genus','genus_taxid','species','species_taxid','scientific_name','taxid','sample']).sum()
    lin_tax_tb.to_csv(f"{tool_path}/{sample}.tsv",sep='\t')
    tables_list.append(lin_tax_tb)

concatDf=pd.concat(tables_list,sort=False,axis=0)
concatDf.to_csv(snakemake.output[0],sep='\t')