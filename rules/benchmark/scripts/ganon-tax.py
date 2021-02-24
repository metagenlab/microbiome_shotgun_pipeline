import pandas as pd
from ete3 import NCBITaxa
import numpy as np
from Bio import Entrez

Entrez.email=snakemake.params.NCBI_email
Entrez.api_key=snakemake.params.NCBI_key
ncbi = NCBITaxa()

def get_lin_tax(tab, ncbi, target_ranks):
    tax = {}
    tab = tab.replace(np.nan, 0)
    taxid_list = tab['taxid']
    percents = tab['read_percent']
    counts = tab['reads_assigned']
    count_dic = dict(zip(taxid_list, counts))
    percent_dic = dict(zip(taxid_list, percents))
    for taxid in taxid_list:
        tax[taxid] = {'taxid': int(taxid), 'read_percent': percent_dic[taxid], 'read_counts': count_dic[taxid]}
        if int(taxid) > 0:
            scientific_name = ncbi.translate_to_names([taxid])[0]
            tax[taxid]['scientific_name'] = scientific_name
            try:
                lineage = ncbi.get_lineage(taxid)
            except ValueError:
                print(f'taxid:{taxid} not found, searching NCBI taxonomy')
                del_entry_name = Entrez.read(Entrez.esummary(db='taxonomy', id=f'{taxid}'))[0]['ScientificName']
                print(f'deleted entry name: {del_entry_name}')
                simple_name = ' '.join(del_entry_name.split(' ')[
                                       0:2])  # Some times, subspecies or a strain number is added to the scientific name, and ete3 cannot find a match
                print(f'first two words of entry name {simple_name}')
                if 'unclassified' in simple_name:
                    print('discarding "unclassified" in name')
                    simple_name = simple_name.split('unclassified ')[1]
                    print(f'final name: {simple_name}')
                if 'incertae' in simple_name:
                    print('discarding "incertae:" in name')
                    simple_name = simple_name.split('incertae')[0].strip()
                    print(f'final name: {simple_name}')  
                updated_taxid = list(ncbi.get_name_translator([simple_name]).values())[0][0]
                lineage = ncbi.get_lineage(updated_taxid)

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
                            continue
                    tax[taxid][rank] = f'{previous_name}_{rank[0:1]}'
        if int(taxid) == 0:
            tax[taxid]['scientific_name'] = 'unclassified'
    df = pd.DataFrame.from_dict(tax, orient='index')
    df = df.replace(np.nan, 'NA')
    return df

target_ranks = ['superkingdom','phylum','order','family','genus','species']


tb_input=pd.read_csv(snakemake.input[0],sep='\t',names=['taxid','reads_assigned','read_percent','reads_assigned_total','reads_uniquely_assigned','rank','name','NA'])
tb_input.replace('unclassified',0,inplace=True)
tb_linear_tax=get_lin_tax(tb_input,ncbi,target_ranks)
tb_linear_tax = tb_linear_tax.groupby(
    ['superkingdom', 'superkingdom_taxid', 'phylum', 'phylum_taxid', 'order', 'order_taxid', 'family', 'family_taxid',
     'genus', 'genus_taxid', 'species', 'species_taxid', 'scientific_name', 'taxid']).sum()
tb_linear_tax['sample']=[snakemake.wildcards.sample]*len(tb_linear_tax)
tb_linear_tax.to_csv(snakemake.output[0], sep='\t')




