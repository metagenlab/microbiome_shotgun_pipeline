import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()

input_file_list=snakemake.input
dic={}
for path in input_file_list:
    sample=path.split('/')[1]
    dic[sample]=path


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

target_ranks = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']
tables_list=[]

for sample in sample_names:
    pathseq_tb=pd.read_csv(dic[sample],sep='\t')
    pathseq_tb.columns = ['taxid', 'taxonomy', 'rank', 'name', 'kingdom', 'score', 'score_normalized', 'reads',
                          'reads_assigned', 'refence_length']
    lin_tax_tb = get_linear_tax(ncbi, target_ranks,  pathseq_tb)
    lin_tax_tb.rename(columns={'read_counts': f'{sample}'})
    tables_list.append(lin_tax_tb)

concatDf=pd.concat(tables_list,sort=False,axis=0,join='outer')
concatDf.to_csv(snakemake.output[0],sep='\t')