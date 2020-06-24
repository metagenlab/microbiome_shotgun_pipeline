import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
superkingdom=snakemake.params.superkingdom
rank=snakemake.params.rank

gold_standard = pd.read_csv(snakemake.input.gold_standard,sep='\t')
superkingdom_gs=gold_standard[gold_standard.superkingdom==superkingdom]




def groupby_samples(table,tool_name,rank):
    if tool_name=='ezvir' or tool_name=='metaphlan':
        sample_tb=table.groupby(['sample',f'{rank}'],as_index=False)['read_counts'].sum()
        sample_tb['read_counts']=sample_tb['read_counts']//2
    else:
        sample_tb = table.groupby(['sample', f'{rank}'], as_index=False)['read_counts'].sum()
    sample_tb['tool'] = [tool_name] * len(sample_tb)
    return sample_tb

subset_gs=groupby_samples(superkingdom_gs,'gold_standard',rank)

files=snakemake.input.tool_out
tables=[]
tool_list=[]
for file in files:
    tool = os.path.split(file)[1].split('.tsv')[0]
    tool_list.append(tool)
    tb=pd.read_csv(file,sep='\t')
    tb=tb[tb.superkingdom==superkingdom]
    tb=groupby_samples(tb,tool,rank)
    tables.append(tb)
all_tools=pd.concat(tables,sort=False)

gs_tool_tb=pd.concat([subset_gs,all_tools],sort=False)


gs_tool_tb.set_index(f'{rank}',inplace=True)
true_genomes_list=list(set(gs_tool_tb.loc[gs_tool_tb['tool']=='gold_standard'].index))
gs_tool_tb.reset_index(inplace=True)
true_pos=gs_tool_tb.loc[gs_tool_tb[f'{rank}'].isin(true_genomes_list)]
true_pos=true_pos.sort_values(by='read_counts',ascending=False)

tool_list.insert(0,'gold_standard')
plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(left=0.4,right=0.81)
bp=sns.barplot(data=true_pos,y=f'{rank}',x='read_counts',hue='tool',hue_order=tool_list)
bp.legend(loc=2, bbox_to_anchor=(1.05, 1),borderaxespad=0.1)
bp.get_figure().savefig(snakemake.output[0])