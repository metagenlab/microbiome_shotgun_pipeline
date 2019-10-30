import pandas as pd

def concat_tables(file_list):
    table_list = []
    for file in file_list:
        table=pd.read_csv(file,sep='\t',index_col=0)
        table.index.name = 'tool'
        table_list.append(table)
    all_tools = pd.concat(table_list, sort=False, axis=0)
    if 'F1_score' in all_tools.columns:
        all_tools=all_tools.sort_values(by='F1_score',ascending=False)
    return all_tools

vir_scores=concat_tables(snakemake.input.vir_scores)
vir_scores.to_csv(snakemake.output.all_vir_scores,sep='\t')
vir_presence=concat_tables(snakemake.input.vir_presence)
vir_presence.to_csv(snakemake.output.all_vir_presence,sep='\t')
bac_scores=concat_tables(snakemake.input.bac_scores)
bac_scores.to_csv(snakemake.output.all_bac_scores,sep='\t')
bac_presence=concat_tables(snakemake.input.bac_presence)
bac_presence.to_csv(snakemake.output.all_bac_presence,sep='\t')

