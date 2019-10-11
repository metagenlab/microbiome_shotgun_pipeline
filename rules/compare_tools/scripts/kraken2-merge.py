import pandas as pd
input_file_list=snakemake.input
sample_names=[name.split("/")[1] for name in input_file_list]
dic=dict(zip(sample_names,input_file_list))

tables_list=[]
for i in sample_names:
    tb=pd.read_csv(dic[i],sep='\t',header=None)
    tb.columns=["percent_covered","rooted_frag_num",f"{i}","rank_code","taxid","scientific_names"]
    sample_df=tb[[f"{i}","taxid"]]
    sample_df.set_index('taxid',drop=True,inplace=True)
    tables_list.append(sample_df)

concatDf=pd.concat(tables_list,axis=1)
concatDf.to_csv(snakemake.output[0],sep='\t')