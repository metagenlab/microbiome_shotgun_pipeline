import pandas as pd
table_list=[]
for file in snakemake.input:
    table=pd.read_csv(file,sep='\t',index_col=0)
    table.index.name = 'superkingdom'
    table_list.append(table)

all_tools=pd.concat(table_list,sort=False,axis=0)
all_tools=all_tools.sort_values(by='F1_score',ascending=False)
all_tools.to_csv(snakemake.output[0],sep='\t')