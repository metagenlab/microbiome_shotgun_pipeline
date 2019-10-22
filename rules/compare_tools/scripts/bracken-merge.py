import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()

input_file_list=snakemake.input
dic={}
for path in input_file_list:
    sample=path.split('/')[1]
    dic[sample]=path

def get_lin_tax(tab,ncbi,target_ranks):
    tax={}
    for taxid in tab['taxid']:
        tax[taxid]={'taxid':taxid,'read_percent':tab.loc[taxid]['read_percent'],'read_counts':tab.loc[taxid]['reads_assigned']}
        try:
            int(taxid)
            lin_txid = ncbi.get_lineage(taxid)
            lin_translation = ncbi.get_taxid_translator(lin_txid)
            taxid2rank = ncbi.get_rank(list(lin_translation.keys()))
            taxid_path=[]
            for n,k in enumerate(taxid2rank.keys()):
                if taxid2rank[k] in target_ranks:
                    tax[taxid][taxid2rank[k]] = lin_translation[k]
                    tax[taxid]['taxid_path']=taxid_path.append(str(lin_txid[n]))
            tax[taxid]['taxid_path']='|'.join(taxid_path)
        except ValueError:
            for j in target_ranks:
                tax[taxid][j]='NA'
                tax[taxid]['taxid_path']='NA'
    df=pd.DataFrame.from_dict(tax,orient='index')
    return df

target_ranks = ['superkingdom','phylum','order','family','genus','species']
tables_list=[]

for sample in dic.keys():
    tb=pd.read_csv(dic[sample],sep='\t')
    tb.columns=['name','taxid','taxonomy_lvl','kraken_assigned_reads','added_reads','reads_assigned','read_percent']
    lin_tax_tb=get_lin_tax(tb,ncbi,target_ranks)
    lin_tax_tb=lin_tax_tb.rename(columns={'read_counts': f'{sample}_counts','read_percent': f'{sample}_percent'})

    lin_tax_tb=lin_tax_tb.groupby(['taxid','species', 'genus', 'family', 'order', 'phylum', 'superkingdom','taxid_path'],as_index=False).sum()
    tables_list.append(lin_tax_tb)

concatDf=pd.concat(tables_list,sort=False,axis=0)
concatDf.to_csv(snakemake.output[0],sep='\t')