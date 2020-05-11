import pandas as pd
import os

tb=pd.read_csv(snakemake.input[0],sep='\t',names=['fasta_header','nucleotide_pos','nucleotide_cov'])
gb=tb.groupby(['fasta_header','nucleotide_pos'],as_index=False)['nucleotide_cov'].sum()
gb.index=gb['fasta_header']
indexes=list(set(gb.index))#list of unique indexes/fasta files
dataframe_dic={}

os.mkdir(snakemake.params[0])
for index in indexes:
    table=gb.loc[gb['fasta_header']==index]
    head=index.split('|')
    acc=head[len(head)-2].replace('.','_')
    path=os.path.join(snakemake.params[0],f"{acc}.csv")
    table.to_csv(path,header=None,index=False,sep='\t')



