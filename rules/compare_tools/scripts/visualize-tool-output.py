import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def concat_sample_counts(file, rank, tool_name, superkingdom, threshold):
    tb_list = []
    col_sel = []
    tool_tb = pd.read_csv(file, sep='\t')
    tool_tb = tool_tb[tool_tb.superkingdom == superkingdom]
    tool_tb = tool_tb.groupby(f'{rank}').sum()

    for col in tool_tb:
        if 'counts' in col:
            col_sel.append(col)
    subset = tool_tb[col_sel]

    for col in subset.columns:
        tb = subset[[col]]
        tb = tb[tb > threshold].dropna()
        type_list = [f'{tool_name}'] * len(tb)
        tb.insert(loc=1, column='tool', value=type_list)
        tb = tb.rename(columns={f"{col}": "read_counts"})
        tb_list.append(tb)

    all_samples = pd.concat(tb_list, axis=0, sort=False)
    all_samples = all_samples.replace(np.nan, 0)
    all_samples.index.name = f'{rank}'
    return all_samples

tool2file={}
outputs=snakemake.input
rank=snakemake.params.rank
value=snakemake.params.value
threshold=snakemake.params.threshold
superkingdom=snakemake.params.superkingdom
for file in outputs:
    tool=file.split('/')[1]
    tool2file[tool]=file
tables=[]
for tool in tool2file.keys():
    tb=concat_sample_counts(tool2file[tool],rank,tool,superkingdom,threshold)
    tables.append(tb)
all_tools = pd.concat(tables, axis=0, sort=False)
all_tools=all_tools.sort_values(by='read_counts',ascending=False)
all_tools=all_tools.reset_index()

plt.figure(figsize=(11.7,8.27))
bp=sns.catplot(data=all_tools,x='read_counts',y=f'{rank}',col='tool',col_wrap=3,kind='bar')
for ax in bp.axes.flat:
    for ylab in ax.get_yticklabels():
        ylab.set_size(6)
    for xlab in ax.get_xticklabels():
        xlab.set_size(6)
        xlab.set_rotation(45)
bp.savefig(snakemake.output.barplot)
pivot=all_tools.pivot_table(all_tools,index=f'{rank}',columns='tool')
pivot=pivot.replace(np.nan,0)
pivot.columns=pivot.columns.droplevel()
pivot=pivot.sort_values(by=list(pivot.columns),ascending=False)


plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(left=0.4, bottom=0.4)
hmc = sns.heatmap(
    pivot,
    square=False,
    cbar_kws={'fraction' : 0.01},
    cmap='OrRd',
    linewidth=1
)
hmc.get_figure().savefig(snakemake.output.heatmap_counts)

def get_presence(table,rank,tool):
    all_found=list(set(table[f'{rank}']))
    tool_found=list(set(table[table['tool']==f'{tool}'][f'{rank}']))
    dic={}
    for hit in all_found:
        dic[hit]={}
        if hit in tool_found:
            dic[hit][f'{tool}']=True
        else:
            dic[hit][f'{tool}']=False
    df=pd.DataFrame.from_dict(dic, orient='index')
    df.index.name=f'{rank}'
    return df

tool_list=list(set(all_tools['tool']))
tb_list=[]
for tool in tool_list:
    tb_list.append(get_presence(all_tools,rank,tool))
matrix=pd.concat(tb_list,axis=1,sort=False)
matrix=matrix.sort_values(by=list(matrix.columns),ascending=False)

plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(left=0.4, bottom=0.4)
categories=sorted(np.unique(matrix))
if len(categories)==1 and categories[0]==True:
    myColors = [(0.3, 0.6, 0.1, 1.0),(0.0, 0.0, 0.0, 0.9)]
else:
    myColors = [(0.0, 0.0, 0.0, 0.9), (0.3, 0.6, 0.1, 1.0)]#RGB color palette
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
hmb = sns.heatmap(
    matrix,
    square=False,
    cbar_kws={'fraction' : 0.01},
    cmap=cmap,
    linecolor='black',linewidths='0.1'
)
colorbar = hmb.collections[0].colorbar
colorbar.set_ticks([0,1.0])
colorbar.set_ticklabels(['absent', 'present'])
hmb.get_figure().savefig(snakemake.output.heatmap_presence)