import pandas as pd

scores_files=snakemake.input.scores
presence_files=snakemake.input.presence


def concat_tables(files):
    tb_list=[]
    for file in files:
        tb=pd.read_csv(file,sep='\t')
        tb_list.append(tb)
    all_tables=pd.concat(tb_list,sort=False,axis=0)
    return all_tables


scores=concat_tables(scores_files)
scores=scores.sort_values(by='F1_score',ascending=False)
scores.to_csv(snakemake.output.all_scores,sep='\t',index=None)
presence=concat_tables(presence_files)
presence.to_csv(snakemake.output.all_presence,sep='\t',index=None)

