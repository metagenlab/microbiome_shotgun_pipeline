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

scores=concat_tables(snakemake.input.scores)
scores.to_csv(snakemake.output.all_scores,sep='\t')
presence=concat_tables(snakemake.input.presence)
presence.to_csv(snakemake.output.all_presence,sep='\t')

