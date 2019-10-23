import pandas as pd
import numpy as np
import re
import math
import copy
from Bio import Entrez
Entrez.email=snakemake.config["email"]
Entrez.api_key=snakemake.config["api_key"]
input_file_list=snakemake.input
sample_names=[name.split("/")[1] for name in input_file_list]
dic=dict(zip(sample_names,input_file_list))

read_length=snakemake.params["read_len"]
target_ranks = ['superkingdom','phylum','order','family','genus','species']

def get_filtered_read_counts(table,read_len,threshold):
    read_counts = {}
    for genome_name in table.index:
        read_num = math.floor(table.loc[genome_name]['TNT'] / read_len)
        if read_num >= threshold:
            read_counts[genome_name] = read_num #dictionary with accessions id as keys and read counts as values
    return read_counts

def get_taxid_from_accession(acc_read_counts_dic):
    acc_id = list(acc_read_counts_dic.keys())
    ids = [re.split(r'_[0-9+]$', acc)[0] for acc in acc_id]#For the EZvir rules, files were named after accession IDs with underscores replacing points. To search for the accession in NCBI we split the underscore and the version number at the end
    uid_list = Entrez.read(Entrez.esearch(db='nucleotide', term=','.join(ids), retmax=100000))['IdList']
    nuc_tax_link = Entrez.read(Entrez.elink(dbfrom='nucleotide', db='taxonomy', id=uid_list))
    acc_id_tax_dic = {}
    for j, name in enumerate(acc_id):
        taxid = nuc_tax_link[j]['LinkSetDb'][0]['Link'][0]['Id']
        acc_id_tax_dic[name] = taxid #dictionary with accessions as keys and taxid as values. Its purpose is to link accessions with read counts to their taxid and taxonomy
    return acc_id_tax_dic

def get_taxonomy_from_taxid(target_ranks,acc_taxid_dic):
    all_taxid = acc_taxid_dic.values()
    taxonomy_summary = Entrez.read(Entrez.efetch(db='taxonomy', id=','.join(all_taxid)))
    tax = {}
    for i in range(len(taxonomy_summary)):
        taxid = taxonomy_summary[i]['TaxId']
        tax[taxid] = {}
        tax[taxid]['ScientificName'] = taxonomy_summary[i]['ScientificName']
        taxid_path = []
        for j in taxonomy_summary[i]['LineageEx']:
            rank = j['Rank']
            if rank in target_ranks:
                names = j['ScientificName']
                tax[taxid][rank] = names
                taxid_path.append(j['TaxId'])
                tax[taxid]['taxid_path'] = '|'.join(taxid_path)
            else:
                tax[taxid][rank] = 'NA'
                tax[taxid]['taxid_path']='NA'
    return tax

def get_tax_table(table,read_len,threshold,target_ranks):
    read_counts_dic=get_filtered_read_counts(table,read_len,threshold)
    acc_taxid_dic=get_taxid_from_accession(read_counts_dic)
    tax=get_taxonomy_from_taxid(target_ranks,acc_taxid_dic)
    dico = {}
    for genome in acc_taxid_dic.keys():
        tax_id = acc_taxid_dic[genome]
        dico[genome] = copy.deepcopy(tax[tax_id])
        dico[genome]['read_counts'] = read_counts_dic[genome]
        dico[genome]['taxid'] = tax_id
    tb = pd.DataFrame.from_dict(dico, orient='index')
    tb=tb.replace(np.nan,'NA')
    return tb

tables_list=[]
for i in sample_names:
    ezvir_tb=pd.read_csv(dic[i],sep=',',index_col='GN')
    ezvir_parsed_tb=get_tax_table(ezvir_tb,read_length,1,target_ranks)#Threshold =1, get all hits with at least 1 read
    ezvir_parsed_tb=ezvir_parsed_tb.rename(columns={'read_counts': f'{i}_counts'})
    samples_ezvir_tb=ezvir_parsed_tb.groupby(['taxid','species','genus','family','order','phylum','superkingdom','taxid_path']).sum()
    samples_ezvir_tb[f'{i}_percent'] = samples_ezvir_tb[f'{i}_counts'] / samples_ezvir_tb[f'{i}_counts'].sum()
    tables_list.append(samples_ezvir_tb)

all_ezvir_tab = pd.concat(tables_list, sort=False, axis=1)
all_ezvir_tab.to_csv(snakemake.output[0],sep='\t')