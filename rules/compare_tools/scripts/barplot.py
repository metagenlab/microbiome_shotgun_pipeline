import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
tool_name=snakemake.wildcards.tool
rank=snakemake.params.rank


def concat_replicate_counts(file, rank, tool_name):
    tb_list = []
    col_sel = []
    tool_tb = pd.read_csv(file, sep='\t')
    tool_tb=tool_tb.groupby(f'{rank}').sum()
    for col in tool_tb:
        if 'counts' in col:
            col_sel.append(col)
    subset = tool_tb[col_sel]

    for col in subset.columns:
        tb = subset[[col]]
        type_list = [f'{tool_name}'] * len(tb)
        tb.insert(loc=1, column='type', value=type_list)
        tb = tb.rename(columns={f"{col}": "read_counts"})
        tb_list.append(tb)

    all_replicates = pd.concat(tb_list, axis=0)
    all_replicates = all_replicates.replace(np.nan, 0)
    return all_replicates


gold_standard = pd.read_csv(snakemake.input.gold_standard,sep='\t')
gold_standard = gold_standard.groupby(f'{rank}').sum()
gold_standard['type'] = ['gold_standard'] * len(gold_standard)# add a column to specify if the read counts are from a tool or the gold standard
true = gold_standard[['read_counts', 'type']]
tool_out = concat_replicate_counts(snakemake.input.tool_out,rank, tool_name)
cat = pd.concat([true, tool_out], axis=0, sort=False)

if rank=='species':
    cat=cat.drop("Homo sapiens")
elif rank=='genus':
    cat=cat.drop("Homo")

true_genomes_list=list(cat.loc[cat['type']=='gold_standard'].index)
cat=cat.reset_index()
subset=cat.loc[cat[f'{rank}'].isin(true_genomes_list)]
subset=subset.sort_values(by='read_counts',ascending=False)
plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(bottom=0.4,right=0.81)
palette={'gold_standard':'C0',f'{tool_name}':'C1'}
sns.set(font_scale=0.9)
bp=sns.barplot(data=subset,x=f'{rank}',y='read_counts',hue='type',palette=palette,order=sorted(true_genomes_list))
bp.legend(loc=2, bbox_to_anchor=(1.05, 1),borderaxespad=0.1)

for tick in bp.xaxis.get_major_ticks():
    tick.label.set_rotation(90)

bp.get_figure().savefig(snakemake.output[0])