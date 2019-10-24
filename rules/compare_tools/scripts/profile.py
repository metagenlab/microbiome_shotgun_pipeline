from ete3 import NCBITaxa
from datetime import date
import pandas as pd
#import numpy as np
today=date.today()
d=today.strftime("%d/%m/%Y")
ncbi = NCBITaxa()
sample=snakemake.wildcards.sample

def get_lin_tax(table,ncbi,target_ranks):
    taxid_list=table['Taxid']
    counts=table['SimReads']
    count_dic = dict(zip(taxid_list, counts))
    tax={}
    for taxid in taxid_list:
        scientific_name=ncbi.translate_to_names([taxid])[0]
        tax[taxid]={'scientific_name':scientific_name,'read_counts':count_dic[taxid]}
        lineage=ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks=ncbi.get_rank(lineage)
        path=[]
        for sub_taxid in lineage:
            rank=ranks[sub_taxid]
            name=names[sub_taxid]
            if rank in target_ranks:
                tax[taxid][rank]=name
                path.append(str(sub_taxid))
        last_taxid=path[len(path)-1]
        last_rank=ranks[int(last_taxid)]
        tax[taxid]['last_taxid']=last_taxid
        tax[taxid]['last_rank']=last_rank
        tax[taxid]['taxid_path']='|'.join(path)
    df=pd.DataFrame.from_dict(tax, orient='index')
    df.index.name='taxid'
    return df
target_ranks = ['superkingdom','phylum','order','family','genus','species']
table=pd.read_csv(snakemake.input[0],sep='\t')
lin_tax_tab=get_lin_tax(table,ncbi,target_ranks)
tb=lin_tax_tab[['last_taxid','last_rank','taxid_path','read_counts']]
tb.columns=['TAXID','RANK','TAXPATH','COUNTS']
tb_human_filtered=tb.drop([9606])#drop human reads
tb_human_filtered['PERCENTAGE']=tb_human_filtered['COUNTS'].div(sum(tb_human_filtered['COUNTS'])).mul(100)
tb_human_filtered=tb_human_filtered[['TAXID','RANK','TAXPATH','PERCENTAGE']]
str_tb=tb_human_filtered.to_string(index=False)

profile_file=open(snakemake.output[0],"w+")
profile_file.write(f'@SampleID:{sample}\n\
@Version:0.9.1"\n\
@Ranks:superkingdom|phylum|order|family|genus|species\n\
@TaxonomyID:ncbi-taxonomy_{d}\n\
@@ {str_tb}')
profile_file.close()