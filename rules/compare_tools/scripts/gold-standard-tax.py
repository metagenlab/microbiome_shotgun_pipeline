from ete3 import NCBITaxa
import pandas as pd
import numpy as np
ncbi = NCBITaxa()

community=snakemake.params.community
def get_lin_tax(table,ncbi,target_ranks,community):
    if 'sample' in table.columns: #Check if input table has a column with sample names in it
        samples=table['sample']
    else: #If the input table does not have a sample column, the sample column is filled with the community name (in config.yaml)
          #For my simulated reads, a metagenome profile was established for each sample.
          #All samples have the same content and the sample column is filled with the name of the meatagenome profile
        samples=[community]*len(table)

    sample_dic = dict(zip(table.index, samples))
    tax={}
    for i in range(len(table)):
        taxid=int(table.iloc[i]['taxid'])
        counts=table.iloc[i]['read_counts']
        scientific_name=ncbi.translate_to_names([taxid])[0]

        tax[i]={'taxid':str(taxid),'scientific_name':scientific_name,'read_counts':counts,'sample':sample_dic[i]}

        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        rank_list = list(ranks.values())
        name_list = list(names.values())
        rank2names = dict(zip(rank_list, name_list))

        for sub_taxid in lineage:
            rank = ranks[sub_taxid]
            if rank in target_ranks:
                tax[i][f"{rank}_taxid"] = sub_taxid
        for n, rank in enumerate(target_ranks):
            if rank in rank2names.keys():
                tax[i][rank] = rank2names[rank]
            else:
                for ranks_to_skip in range(1, len(target_ranks) + 1):
                    previous_rank = target_ranks[n - ranks_to_skip]
                    if previous_rank in rank2names.keys():
                        previous_name = rank2names[previous_rank]
                        break
                    if previous_rank not in rank2names.keys():
                        continue
                tax[i][rank] = f'{previous_name}_{rank[0:1]}'
    df = pd.DataFrame.from_dict(tax, orient='index')
    df = df.replace(np.nan, 'NA')
    return df
target_ranks = ['superkingdom','phylum','order','family','genus','species']
table=pd.read_csv(snakemake.input[0],sep='\t')
rename_dic={'Taxid':'taxid','SimReads':'read_counts'}
table.rename(rename_dic,axis='columns',inplace=True)

lin_tax_tab=get_lin_tax(table,ncbi,target_ranks,community)
column_order=['superkingdom','superkingdom_taxid','phylum','phylum_taxid','order','order_taxid','family','family_taxid','genus','genus_taxid','species','species_taxid','scientific_name','taxid','sample','read_counts']
lin_tax_tab=lin_tax_tab[column_order]
lin_tax_tab.to_csv(snakemake.output[0],sep='\t',index=None)
