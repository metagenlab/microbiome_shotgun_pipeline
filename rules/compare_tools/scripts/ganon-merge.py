import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()



virus_input=snakemake.input['vir']
bacteria_input=snakemake.input['bac']

vir_dic={}
for path in virus_input:
    sample=path.split('/')[1]
    vir_dic[sample]=path

bac_dic={}
for path in bacteria_input:
    sample=path.split('/')[1]
    bac_dic[sample]=path


def get_linear_tax(ncbi,target_ranks,table):
    taxid_read_counts_dic = dict(zip(table['taxid'], table['reads_assigned']))
    taxid_names_dic = dict(zip(table['taxid'], table['name']))
    taxid_list=list(taxid_read_counts_dic.keys())
    tax = {}
    for taxid in taxid_list:
        tax[taxid] = {'taxid': taxid, 'name': taxid_names_dic[taxid], 'read_counts': taxid_read_counts_dic[taxid]}
        try:
            int(taxid)
            lin_txid = ncbi.get_lineage(taxid)
            lin_translation = ncbi.get_taxid_translator(lin_txid)
            taxid2rank = ncbi.get_rank(list(lin_translation.keys()))
            for k in taxid2rank.keys():
                if taxid2rank[k] in target_ranks:
                    tax[taxid][taxid2rank[k]] = lin_translation[k]
        except ValueError:#If the taxid is not an integer
            for j in target_ranks:
                tax[taxid][j] = 'NA'
    linear_tax_tab=pd.DataFrame.from_dict(tax,orient='index')
    return linear_tax_tab

sample_table_list=[]
target_ranks = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']

for sample in vir_dic.keys():
    vir_tb_input=pd.read_csv(vir_dic[sample],sep='\t',names=['taxid','reads_assigned','read_percent','reads_assigned_total','reads_uniquely_assigned','rank','name','NA'])
    vir_tb_linear_tax=get_linear_tax(ncbi,target_ranks,vir_tb_input)
    bac_tb_input=pd.read_csv(bac_dic[sample],sep='\t',names=['taxid','reads_assigned','read_percent','reads_assigned_total','reads_uniquely_assigned','rank','name','NA'])
    bac_tb_linear_tax=get_linear_tax(ncbi,target_ranks,bac_tb_input)
    all_tb=pd.concat([bac_tb_linear_tax, vir_tb_linear_tax], sort=False)
    all_grouped=all_tb.groupby(['taxid', 'name', 'species', 'genus', 'family', 'order', 'phylum', 'superkingdom']).sum()
    all_grouped.rename(columns={'read_counts':f'{sample}'})
    sample_table_list.append(all_grouped)

all_ganon=pd.concat(sample_table_list,sort=False,axis=0,join='outer')
all_ganon.to_csv(snakemake.output[0],sep='\t')