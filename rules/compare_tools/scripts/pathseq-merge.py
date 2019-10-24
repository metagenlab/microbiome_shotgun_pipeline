import pandas as pd
from ete3 import NCBITaxa
import numpy as np
ncbi = NCBITaxa()

input_file_list=snakemake.input
dic={}
for path in input_file_list:
    sample=path.split('/')[1]
    dic[sample]=path


def get_lin_tax(tb, ncbi, target_ranks):

    tb = tb.loc[tb['reads_assigned'] > 0]
    tax_path = tb['taxonomy']
    counts = tb['reads_assigned']
    taxid_list = tb['taxid']
    count_dic = dict(zip(taxid_list, counts))
    path_dic = dict(zip(taxid_list, tax_path))

    tax = {}
    for taxid in list(count_dic.keys()):
        tax[taxid] = {'read_counts': count_dic[taxid]}
        if int(taxid) > 1:
            try:
                tax[taxid] = {'taxid': taxid, 'read_counts': count_dic[taxid]}
                lineage = ncbi.get_lineage(taxid)
            except ValueError:
                taxonomy_path = path_dic[taxid]
                split_path = taxonomy_path.split('|')
                name_rank_before = split_path[len(split_path) - 2]
                newtaxid = list(ncbi.get_name_translator([name_rank_before]).values())[0][0]
                tax[taxid] = {'taxid': newtaxid, 'read_counts': count_dic[taxid]}
                lineage = ncbi.get_lineage(newtaxid)

            names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)
            scientific_name = ncbi.translate_to_names([taxid])[0]
            tax[taxid]['scientific_name'] = scientific_name

            path = []
            for sub_taxid in lineage:
                rank = ranks[sub_taxid]
                name = names[sub_taxid]
                if rank in target_ranks:
                    tax[taxid][rank] = name
                    path.append(str(sub_taxid))

            if len(path) > 0:
                last_taxid = path[len(path) - 1]
                last_rank = ranks[int(last_taxid)]
                tax[taxid]['last_taxid'] = last_taxid
                tax[taxid]['last_rank'] = last_rank
                tax[taxid]['taxid_path'] = '|'.join(path)
        else:
            tax[taxid]['scientific_name'] = np.nan
            tax[taxid]['last_taxid'] = np.nan
            tax[taxid]['last_rank'] = np.nan
            tax[taxid]['taxid_path'] = np.nan

    df = pd.DataFrame.from_dict(tax, orient='index')
    return df
target_ranks = ['superkingdom','phylum','order','family','genus','species']
tables_list=[]

for sample in dic.keys():
    pathseq_tb=pd.read_csv(dic[sample],sep='\t')
    pathseq_tb.columns = ['taxid', 'taxonomy', 'rank', 'name', 'kingdom', 'score', 'score_normalized', 'reads','reads_assigned', 'refence_length']
    lin_tax_tb = get_lin_tax(pathseq_tb,ncbi,target_ranks)

    lin_tax_tb=lin_tax_tb.rename(columns={'read_counts': f'{sample}_counts'})

    lin_tax_tb[f'{sample}_percent'] = lin_tax_tb[f'{sample}_counts'].div(sum(lin_tax_tb[f'{sample}_counts'])).mul(100)

    lin_tax_tb = lin_tax_tb.groupby(['taxid', 'scientific_name', 'species', 'genus', 'family', 'order', 'phylum', 'superkingdom', 'taxid_path','last_taxid', 'last_rank']).sum()
    lin_tax_tb.to_csv(f"benchmark_tools/tables/pathseq/{sample}.tsv", sep='\t')
    tables_list.append(lin_tax_tb)

concatDf=pd.concat(tables_list,sort=False,axis=1)
concatDf.to_csv(snakemake.output[0],sep='\t')