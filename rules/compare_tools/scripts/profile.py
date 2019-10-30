from ete3 import NCBITaxa
import pandas as pd
import numpy as np
ncbi = NCBITaxa()
def get_lin_tax(table,ncbi,target_ranks):
    taxid_list=table['taxid']
    counts=table['read_counts']
    count_dic = dict(zip(taxid_list, counts))
    tax={}
    for taxid in taxid_list:
        scientific_name=ncbi.translate_to_names([taxid])[0]
        tax[taxid]={'taxid':str(taxid),'scientific_name':scientific_name,'read_counts':count_dic[taxid]}
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
                for ranks_to_skip in range(1, len(target_ranks) + 1):
                    previous_rank = target_ranks[n - ranks_to_skip]
                    if previous_rank in rank2names.keys():
                        previous_name = rank2names[previous_rank]
                        break
                    if previous_rank not in rank2names.keys():
                        continue
                tax[taxid][rank] = f'{previous_name}_{rank[0:1]}'
    df = pd.DataFrame.from_dict(tax, orient='index')
    df = df.replace(np.nan, 'NA')
    return df
target_ranks = ['superkingdom','phylum','order','family','genus','species']
table=pd.read_csv(snakemake.input[0],sep='\t')
table.columns=['names','read_counts','taxid','superkingdom','phylum','order','family','genus','species']
lin_tax_tab=get_lin_tax(table,ncbi,target_ranks)
lin_tax_tab.to_csv(snakemake.output[0],sep='\t',index=None)
