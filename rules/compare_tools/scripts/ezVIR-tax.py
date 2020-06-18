import pandas as pd
import math
import copy
from ete3 import NCBITaxa
import numpy as np

ncbi = NCBITaxa()

input_file=snakemake.input.hits
genome_taxids_path=snakemake.params.gen_taxids
read_length=snakemake.params.read_len
target_ranks = ['superkingdom','phylum','order','family','genus','species']


def get_filtered_read_counts(table,threshold):
    read_counts = {}
    for genome_name in table.index:
        read_num = table.loc[genome_name]['TMR']
        if read_num >= threshold:
            read_counts[genome_name] = read_num #dictionary with accessions id as keys and read counts as values
    return read_counts

def get_taxid_from_accession(acc_read_counts_dic,taxid_tb):
    acc_id = list(acc_read_counts_dic.keys())
    acc_id_tax_dic = {}
    for acc in acc_id:
        taxid = taxid_tb.loc[acc]['taxid']
        acc_id_tax_dic[acc] = taxid #dictionary with accessions as keys and taxid as values. Its purpose is to link accessions with read counts to their taxid and taxonomy
    return acc_id_tax_dic

def get_taxonomy_from_taxid(target_ranks,acc_taxid_dic):
    tax={}
    taxid_list = list(acc_taxid_dic.values())
    for taxid in taxid_list:
        tax[taxid] = {'taxid':str(taxid)}
        if int(taxid) > 0:

            scientific_name = list(ncbi.get_taxid_translator([taxid]).values())
            if len(scientific_name) == 0:
                tax[taxid]['scientific_name'] = 'NA'
            else:
                tax[taxid]['scientific_name'] = str(scientific_name[0])

            lineage = ncbi.get_lineage(taxid)
            names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)
            rank_list = list(ranks.values())
            name_list = list(names.values())
            rank2names = dict(zip(rank_list, name_list))

            ranks_found=[]
            for sub_taxid in lineage:
                rank = ranks[sub_taxid]
                ranks_found.append(rank)
                if rank in target_ranks:
                    tax[taxid][f"{rank}_taxid"] = sub_taxid
            missing_ranks=[x for x in target_ranks if x not in ranks_found]
            for m_rank in missing_ranks:
                tax[taxid][f"{m_rank}_taxid"]= 'NA'


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
    return tax

def get_tax_table(table,taxid_tb,threshold,target_ranks):
    read_counts_dic=get_filtered_read_counts(table,threshold)
    acc_taxid_dic=get_taxid_from_accession(read_counts_dic,taxid_tb)
    tax=get_taxonomy_from_taxid(target_ranks,acc_taxid_dic)
    dico = {}
    for genome in acc_taxid_dic.keys():
        tax_id = acc_taxid_dic[genome]
        dico[genome] = copy.deepcopy(tax[tax_id])
        dico[genome]['read_counts'] = read_counts_dic[genome]
    tb = pd.DataFrame.from_dict(dico, orient='index')
    tb=tb.replace(np.nan,'NA')
    return tb


genome_taxids_tb=pd.read_csv(genome_taxids_path,sep=',')
genome_taxids_tb.columns=['acc','taxid']
genome_taxids_tb.set_index('acc',inplace=True)

ezvir_tb=pd.read_csv(input_file,sep=',',index_col='GN')
ezvir_parsed_tb=get_tax_table(ezvir_tb,genome_taxids_tb,1,target_ranks)#Threshold =1, get all hits with at least 1 read
sample_ezvir_tb = ezvir_parsed_tb.groupby(
    ['superkingdom', 'superkingdom_taxid', 'phylum', 'phylum_taxid', 'order', 'order_taxid', 'family', 'family_taxid',
     'genus', 'genus_taxid', 'species', 'species_taxid', 'scientific_name', 'taxid']).sum()
sample_ezvir_tb['sample']=[snakemake.wildcards.sample]*len(sample_ezvir_tb)
sample_ezvir_tb[f'read_percent'] = sample_ezvir_tb[f'read_counts'].div(sum(sample_ezvir_tb[f'read_counts'])).mul(100)
sample_ezvir_tb.to_csv(snakemake.output[0], sep='\t')
