
import pandas as pd
samples_tab=snakemake.input
def concat_tables(files):
    tb_list=[]
    for file in files:
        tb=pd.read_csv(file,sep='\t')
        tb_list.append(tb)
    all_tables=pd.concat(tb_list,sort=False,axis=0)
    return all_tables

cat=concat_tables(samples_tab)
out=cat.to_csv(snakemake.output[0],sep='\t',index=None)