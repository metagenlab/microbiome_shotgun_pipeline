input_file_list=snakemake.input
sample_names=[name.split("/")[1] for name in input_file_list]
dic=dict(zip(sample_names,input_file_list))
tables_list=[]
for i in sample_names:
    pathseq_tb=pd.read_csv(dic[i],sep='\t',index_col='tax_id')
    tb=pathseq_tb[["unambiguous"]]
    tb.columns=[f'{i}']
    tables_list.append(tb)
concatDf=pd.concat(tables_list,axis=1)
concatDf.to_csv(snakemake.output[0],sep='\t')