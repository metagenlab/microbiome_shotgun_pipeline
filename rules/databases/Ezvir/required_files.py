from Bio import SeqIO
from Bio import Entrez
import pandas as pd
from ete3 import NCBITaxa
import random
from random import randint
import logging
logging.basicConfig(format='%(asctime)s %(message)s',filename=snakemake.log[0], level=logging.DEBUG)


ncbi = NCBITaxa()
update=False
if update == True:
    ncbi.update_taxonomy_database()



Entrez.email = snakemake.params.email
Entrez.api_key = snakemake.params.api_key
input_fna = snakemake.input.fna
input_acc2taxid = snakemake.input.acc2taxid
output_genome_lengths = snakemake.output.len
output_acc2taxid = snakemake.output.taxids
output_genome_names = snakemake.output.names
output_colours = snakemake.output.colours
random.seed(0)  # to fix virus family color array


def get_taxonomy_names(taxid, target_ranks):
    try:
        common_name_dict = ncbi.get_taxid_translator([taxid])
        common_name = list(common_name_dict.values())[0]
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        rank_list = list(ranks.values())
        name_list = list(names.values())
    except ValueError:
        rank_list = 'na'
        name_list = 'na'
        common_name = 'na'
    rank2names = dict(zip(rank_list, name_list))
    dic = {}
    dic['common_name'] = common_name
    for rank in target_ranks:
        if rank in rank2names.keys():
            dic[rank] = rank2names[rank]
        else:
            dic[rank] = 'na'
    return dic


handle = SeqIO.parse(input_fna, 'fasta')
acc_list = []
acc_length = {}
acc_desc = {}
acc_no_version = []
for record in handle:
    desc = record.description.split(record.id + ' ')[1].replace(' ', '_').replace(',', '-').replace('/', '-').replace(
        '[', '_').replace(']', '_').replace('{', '_').replace('}', '_')
    acc_split = record.id.split('|')
    acc_no_version.append(acc_split[len(acc_split) - 2].split('.')[0])
    acc = acc_split[len(acc_split) - 2].replace('.', '_')
    acc_list.append(acc)
    acc_length[acc] = len(record.seq)
    acc_desc[acc] = desc

acc2len_tb = pd.DataFrame.from_dict(acc_length, orient='index', columns=['genome_length'])
acc2len_tb.to_csv(output_genome_lengths, sep='\t')
logging.info('genome lengths: Done')
acc_list_sorted = sorted(acc_list)
acc_list_sorted_no_ver = sorted(acc_no_version)

ver_to_acc = dict(zip(acc_list_sorted, acc_list_sorted_no_ver))#Dictionary that links accessions with versions to acessions without versions

acc2taxid = pd.read_csv(input_acc2taxid, sep='\t', chunksize=1000000)
ver_list = []

for df_chunk in acc2taxid:
    ver_list.append(df_chunk[df_chunk["accession"].isin(acc_list_sorted_no_ver)])

acc_taxid_tb = pd.concat(ver_list)


acc_to_taxid_found = dict(zip(acc_taxid_tb['accession'], acc_taxid_tb['taxid']))

logging.info(f'{len(acc_to_taxid_found)} entries for which there is a taxid')

removed_acc=list(set(acc_no_version)-set(acc_to_taxid_found.keys()))#Number of entries for which no matching taxid was found

logging.info(f'{len(removed_acc)} entries for which no taxid was found')

##Search Entrez for removed nucleotide entries
rm_acc2taxid={}
for rm_acc in removed_acc:
    try:
        taxid=int(Entrez.read(Entrez.esummary(db='nucleotide', id=rm_acc))[0]['TaxId'])
        logging.info(f'taxid {taxid} found for removed entry {rm_acc}')
    except RuntimeError:
        taxid='NA'
        logging.info(f'taxid not found for removed entry {rm_acc}')
    rm_acc2taxid[rm_acc]=taxid

accessions2taxid={}
for acc in ver_to_acc:
    acc_no_ver=ver_to_acc[acc]
    try:
        accessions2taxid[acc]=acc_to_taxid_found[acc_no_ver]#Try to find corresponding taxid in accessions to taxid dictionary
    except KeyError:
        accessions2taxid[acc]=rm_acc2taxid[acc_no_ver]#Try to find corresponding taxid in removed accessions to taxid dictionary


acc2taxid_tb=pd.DataFrame.from_dict(accessions2taxid,orient='index',columns=['taxid'])
acc2taxid_tb.index.name='accession'
acc2taxid_tb.to_csv(output_acc2taxid,sep=',')
logging.info('genome taxids: Done')

unique_taxid_list = list(set(acc2taxid_tb['taxid']))

taxid2ranknames={}#dictionary listing viral families for each taxid
for taxid in unique_taxid_list:
    info=get_taxonomy_names(taxid,['family'])
    taxid2ranknames[taxid]=info


genome_names={}
acc_list=list(accessions2taxid.keys())
for acc in enumerate(acc_list):
    taxid=accessions2taxid[acc]
    genome_names[acc]={}
    genome_names[acc]['header_name']=acc_desc[acc]
    genome_names[acc]['family']=taxid2ranknames[taxid]['family']
    genome_names[acc]['species']=taxid2ranknames[taxid]['common_name']

genome_names_tb=pd.DataFrame.from_dict(genome_names,orient='index')
genome_names_tb.to_csv(output_genome_names,sep=',',header=None)
logging.info('genome names: Done')

unique_families=list(set(genome_names_tb['family']))
colors=[]
for i in range(len(unique_families)):
    colors.append('#%06X' % randint(0, 0xFFFFFF))

family2color=dict(zip(unique_families,colors))
color_tb=pd.DataFrame.from_dict(family2color,orient='index',columns=['colorname'])
color_tb.index.name='v_grp'
color_tb.to_csv(output_colours,sep=',')
logging.info('genome families colour: Done')