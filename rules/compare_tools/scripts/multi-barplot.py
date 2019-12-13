import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

superkingdom=snakemake.params.superkingdom
rank=snakemake.params.rank

gold_standard = pd.read_csv(snakemake.input.gold_standard,sep='\t')
superkingdom_gs=gold_standard[gold_standard.superkingdom==superkingdom]

superkingdom_gs=superkingdom_gs.groupby(f'{rank}',as_index=False).sum()
superkingdom_gs['tool']=['gold_standard']*len(superkingdom_gs)
subset_gs=superkingdom_gs[[f'{rank}','read_counts','tool']]



def groupby_samples(table,tool_name,rank):
    samples=sorted(list(set(table['sample'])))
    tables=[]

    for sample in samples:
        sample_tb = table[table['sample'] == sample]
        sample_tb = sample_tb.groupby(f'{rank}', as_index=False).sum()
        sample_tb['tool'] = [tool_name] * len(sample_tb)
        subset=sample_tb[[f'{rank}','read_counts','tool']]
        tables.append(subset)
    all_samples=pd.concat(tables,sort=False)
    return all_samples


files=snakemake.input.tool_out
tables=[]
tool_list=[]
for file in files:
    tool = file.split('/')[1]
    tool_list.append(tool)
    tb=pd.read_csv(file,sep='\t')
    tb=tb[tb.superkingdom==superkingdom]
    tb=groupby_samples(tb,tool,rank)
    tables.append(tb)
all_tools=pd.concat(tables,sort=False)

gs_tool_tb=pd.concat([subset_gs,all_tools],sort=False)

gs_tool_tb.set_index(f'{rank}',inplace=True)
true_genomes_list=list(gs_tool_tb.loc[gs_tool_tb['tool']=='gold_standard'].index)
gs_tool_tb.reset_index(inplace=True)
true_pos=gs_tool_tb.loc[gs_tool_tb[f'{rank}'].isin(true_genomes_list)]
true_pos=true_pos.sort_values(by='read_counts',ascending=False)

tool_list.insert(0,'gold_standard')
plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(left=0.4,right=0.81)
bp=sns.barplot(data=true_pos,y=f'{rank}',x='read_counts',hue='tool',hue_order=tool_list)
bp.legend(loc=2, bbox_to_anchor=(1.05, 1),borderaxespad=0.1)
bp.get_figure().savefig(snakemake.output[0])