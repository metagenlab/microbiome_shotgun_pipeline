import pandas as pd
from ete3 import NCBITaxa
import numpy as np
ncbi = NCBITaxa()


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
                tax[taxid] = {'taxid': str(taxid), 'read_counts': count_dic[taxid]}
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
                            continue

                    tax[taxid][rank] = f'{previous_name}_{rank[0:1]}'
    df = pd.DataFrame.from_dict(tax, orient='index')
    df=df.replace(np.nan,'NA')
    return df

target_ranks = ['superkingdom','phylum','order','family','genus','species']

pathseq_tb=pd.read_csv(snakemake.input[0],sep='\t',names=['taxid', 'taxonomy', 'rank', 'name', 'kingdom', 'score', 'score_normalized', 'reads','reads_assigned', 'refence_length'],  skiprows=1)

lin_tax_tb = get_lin_tax(pathseq_tb,ncbi,target_ranks)
lin_tax_tb = lin_tax_tb.groupby(
    ['superkingdom', 'superkingdom_taxid', 'phylum', 'phylum_taxid', 'order', 'order_taxid', 'family', 'family_taxid',
     'genus', 'genus_taxid', 'species', 'species_taxid', 'scientific_name', 'taxid']).sum()
lin_tax_tb['sample']=[snakemake.wildcards.sample]*len(lin_tax_tb)
lin_tax_tb['read_percent'] = lin_tax_tb['read_counts'].div(sum(lin_tax_tb['read_counts'])).mul(100)
lin_tax_tb.to_csv(snakemake.output[0], sep='\t')