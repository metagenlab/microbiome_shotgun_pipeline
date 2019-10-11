input_file_list=snakemake.input
sample_names=[name.split("/")[1] for name in input_file_list]
dic=dict(zip(sample_names,input_file_list))
tables_list=[]
for i in sample_names:
    tb=pd.read_csv(dic[i],sep='\t',index_col='taxon_id')
    tb=tb[["reads"]]
    tb.columns=[f'{i}']
    tables_list.append(tb)
concatDf=pd.concat(tables_list,axis=1)
concatDf.to_csv(snakemake.output[0],sep='\t')