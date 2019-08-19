import pandas as pd
tb_krona=pd.read_csv(snakemake.input[0],sep='\t')
tax_list=list(tb_krona['taxonomy'])
tax_list_rep=[tax.replace('|','\t') for tax in tax_list]
read_list=list(tb_krona['unambiguous'])#Only reads with high scores
krona_out=pd.DataFrame(tax_list_rep,read_list)
krona_out.to_csv(snakemake.output[0],sep='\t',header=None)