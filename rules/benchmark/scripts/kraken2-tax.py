import pandas as pd
from Bio import Entrez
from ete3 import NCBITaxa
import numpy as np
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)
logging.captureWarnings(True)
Entrez.email = snakemake.params.NCBI_email
Entrez.api_key = snakemake.params.NCBI_key
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
        if int(taxid) > 0 and count_dic[taxid] > 0:
            tax[taxid] = {'taxid': str(int(taxid)), 'read_percent': percent_dic[taxid],
                          'read_counts': str(int(count_dic[taxid]))}
            lineage=None
            scientific_name = ncbi.translate_to_names([taxid])[0]
            tax[taxid]['scientific_name'] = scientific_name
            try:
                lineage = ncbi.get_lineage(taxid)
                foundlineage = True
            except ValueError:
                logging.error(f'taxid:{taxid} not found, searching NCBI taxonomy')
                try:
                    del_entry_name = Entrez.read(Entrez.esummary(db='taxonomy', id=f'{taxid}'))[0]['ScientificName']
                    logging.info(f'deleted entry name: {del_entry_name}')
                    simple_name = ' '.join(del_entry_name.split(' ')[
                                           0:2])
                    # Some times, subspecies is added to the scientific name,and ete3 cannot find a match
                    logging.info(f'first two words of entry name {simple_name}')
                    if 'unclassified' in simple_name:
                        logging.info('discarding "unclassified" in name')
                        simple_name = simple_name.split('unclassified ')[1]
                        logging.info(f'final name: {simple_name}')
                    if 'incertae' in simple_name:
                        logging.info('discarding "incertae:" in name')
                        simple_name = simple_name.split('incertae')[0].strip()
                        logging.info(f'final name: {simple_name}')
                    updated_taxid = list(ncbi.get_name_translator([simple_name]).values())[0][0]
                    logging.info(f"updated_taxid: {updated_taxid} ")
                    lineage = ncbi.get_lineage(updated_taxid)
                    foundlineage = True
                except RuntimeError:
                    logging.error(f"Could not convert to ncbi taxid, skipping:{tax[taxid]}")
                    foundlineage = False

            if foundlineage:
                names = ncbi.get_taxid_translator(lineage)
                ranks = ncbi.get_rank(lineage)
                rank_list = list(ranks.values())
                name_list = list(names.values())
                rank2names = dict(zip(rank_list, name_list))
                for sub_taxid in lineage:
                    rank = ranks[sub_taxid]
                    if rank in target_ranks:
                        tax[taxid][f"{rank}_taxid"] = str(int(sub_taxid))
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
        if int(taxid) == 0 and count_dic[taxid] > 0:
            tax[taxid] = {'taxid': str(int(taxid)), 'read_percent': percent_dic[taxid],
                          'read_counts': str(int(count_dic[taxid]))}
            tax[taxid]['scientific_name'] = 'unclassified'
    df = pd.DataFrame.from_dict(tax, orient='index')
    return df


target_ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
sample_tb = pd.read_csv(snakemake.input.tax,sep='\t',names=["read_percent","rooted_frag_num","reads_assigned",
                                                         "rank_code","taxid","scientific_names"])
sample_types = pd.read_csv(snakemake.input.sample_info,sep='\t')[['SampleName','bodysite','type']]
sample_types = sample_types.rename(columns={'SampleName':'sample'})
lin_tax_tb = get_lin_tax(sample_tb,ncbi,target_ranks)
lin_tax_tb = lin_tax_tb.astype({'read_counts':'int64','taxid':'object'})
lin_tax_tb['sample']=[snakemake.wildcards.sample]*len(lin_tax_tb)
addtypes=pd.merge(lin_tax_tb, sample_types, how='left', on='sample')

linear_taxonomy_table = addtypes[['superkingdom','superkingdom_taxid','phylum','phylum_taxid','order',
                                  'order_taxid','family','family_taxid','genus','genus_taxid','species','species_taxid',
                                  'strain', 'strain_taxid','scientific_name','taxid','read_counts','read_percent',
                                  'sample','bodysite','type']]
linear_taxonomy_table.replace(np.nan, 'na', inplace=True)
linear_taxonomy_table.to_csv(snakemake.output[0], sep='\t',index=False)
