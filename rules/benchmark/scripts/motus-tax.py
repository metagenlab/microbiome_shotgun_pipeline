import pandas as pd
from Bio import Entrez
Entrez.email=snakemake.params.NCBI_email
Entrez.api_key=snakemake.params.NCBI_key
from ete3 import NCBITaxa
import numpy as np
ncbi = NCBITaxa()

sample=snakemake.wildcards.sample
tsv=pd.read_csv(snakemake.input[0],sep='\t')
tsv=tsv.replace('-1','unclassified')
tsv=tsv.replace(np.nan,0)
tsv['percent']=tsv[sample]/sum(tsv[sample])
tsv['name']=[n.split(' [')[0] for n in tsv['#consensus_taxonomy']]#clean names, i.e, remove motu id number and brackets
tab=tsv[(tsv[sample]>0)]#only get hits with reads


def get_lin_tax(tab, ncbi, target_ranks, sample):
    tax = {}
    for i in range(len(tab)):
        name = tab.iloc[i]['name']
        taxid = tab.iloc[i]['NCBI_tax_id']
        reads = tab.iloc[i][sample]
        percent = tab.iloc[i]['percent']
        if int(taxid) == 0 and name != 'unclassified':#get hits that are classified, but do not have a taxid from motus table
            try:
                taxid = ncbi.get_name_translator([name])[name][0]#try to search for a taxid from motus names
            except KeyError:
                name = name.split(' ')[0]#take genus name if no taxid was found with original species name
                taxid = ncbi.get_name_translator([name])[name][0]

        tax[taxid] = {'taxid': int(taxid), 'read_percent': percent, 'read_counts': reads}
        if int(taxid) > 0:
            scientific_name = ncbi.translate_to_names([int(taxid)])[0]
            tax[taxid]['scientific_name'] = scientific_name
            try:
                lineage = ncbi.get_lineage(taxid)
            except ValueError:
                print(f'taxid:{taxid} not found, searching NCBI taxonomy')
                del_entry_name = Entrez.read(Entrez.esummary(db='taxonomy', id=f'{taxid}'))[0]['ScientificName']
                print(f'deleted entry name: {del_entry_name}')
                simple_name = ' '.join(del_entry_name.split(' ')[
                                       0:2])  # Some times, subspecies is added to the scientific name, and ete3 cannot find a match
                print(f'first two words of entry name {simple_name}')
                if 'unclassified' in simple_name:
                    print('discarding "unclassified" in name')
                    simple_name = simple_name.split('unclassified ')[1]
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
                    tax[taxid][f"{rank}_taxid"] = str(sub_taxid)
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

        if int(taxid) == 0 and name == 'unclassified':
            tax[taxid]['scientific_name'] = 'unclassified'

    df = pd.DataFrame.from_dict(tax, orient='index')
    df = df.replace(np.nan, 'NA')
    return df

target_ranks = ['superkingdom','phylum','order','family','genus','species']
lin_tax=get_lin_tax(tab,ncbi,target_ranks,sample)
lin_tax=lin_tax.groupby(['superkingdom','superkingdom_taxid','phylum','phylum_taxid','order','order_taxid','family','family_taxid','genus','genus_taxid','species','species_taxid','scientific_name','taxid']).sum()
lin_tax['sample']=[snakemake.wildcards.sample]*len(lin_tax)
lin_tax.to_csv(snakemake.output[0],sep='\t')