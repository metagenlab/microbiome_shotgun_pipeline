import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from numpy import mean
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages


def concat_sample_counts(file, target_val, sample_name, rank, tool_name, superkingdom, threshold):
    tb_list = []
    col_sel = []
    tool_tb_raw = pd.read_csv(file, sep='\t')
    tool_tb_raw = tool_tb_raw[tool_tb_raw.superkingdom.str.lower() == superkingdom]
    tool_tb = tool_tb_raw.groupby(f'{rank}').sum()
    for col in tool_tb:
        if target_val in col:
            col_sel.append(col)
    subset = tool_tb[col_sel]
    for col in subset.columns:
        tb = subset[subset > threshold].dropna()
        type_list = [f'{tool_name}'] * len(tb)
        sample_name_list = [f'{sample_name}'] * len(tb)
        tb.insert(loc=1, column='tool', value=type_list)
        tb.insert(loc=2, column='sample', value=sample_name_list)
        tb = tb.rename(columns={f"{col}": f"read_{target_val}"})
        tb_list.append(tb)

    all_replicates = pd.concat(tb_list, axis=0, sort=False)
    all_replicates = all_replicates.replace(np.nan, 0)
    all_replicates.index.name = f'{rank}'
    return all_replicates

outputs=snakemake.input
rank=snakemake.params.rank
value=snakemake.params.value
superkingdom=snakemake.params.superkingdom
path=snakemake.params.matrix_path
threshold=snakemake.params.threshold
sample_list=[]
tables=[]

for file in outputs:
    tool=file.split('/')[1]
    sample=file.split('/')[2].split('.tsv')[0]
    tb=concat_sample_counts(file,value,sample,rank,tool,superkingdom,threshold)
    tables.append(tb)

all_tools = pd.concat(tables, axis=0, sort=False)
all_tools=all_tools.sort_values(by=f'read_{value}',ascending=False)
all_tools=all_tools.reset_index()


unique_sample_list=list(set(all_tools['sample']))
unique_names=list(set(all_tools[f'{rank}']))
square=False
annot=True
fraction = 0.9  # For the presence heatmap: only consider hits identified by 90% of the tools


sns.set(font_scale=1.0,rc={'figure.figsize':(11.7,8.27)})

bp=sns.catplot(data=all_tools,x=f'read_{value}',y=f'{rank}',col='sample',col_wrap=3,kind='bar',estimator=mean,errwidth=0.3)
bp.set(xscale='log')
for ax in bp.axes.flat:
    for ylab in ax.get_yticklabels():
        ylab.set_size(5)
    for xlab in ax.get_xticklabels():
        xlab.set_size(8)
bp.savefig(snakemake.output.barplot_samples)

bpt=sns.catplot(data=all_tools,x=f'read_{value}',y=f'{rank}',hue='sample',col='tool',col_wrap=3,kind='bar',estimator=mean,errwidth=0.3,legend=False)
plt.legend(bbox_to_anchor=(1.10, 1.0),borderaxespad=0)
bpt.set(xscale='log')
for ax in bpt.axes.flat:
    for ylab in ax.get_yticklabels():
        ylab.set_size(5)
    for xlab in ax.get_xticklabels():
        xlab.set_size(8)
bpt.savefig(snakemake.output.barplot_tools)




def draw_heatmap_counts(data,sample,rank,square,annot,tab_path):
    subset=data[data['sample']==sample]
    pivot=subset.pivot_table(subset,index=f'{rank}',columns='tool')
    pivot=pivot.replace(np.nan,0)
    pivot.columns=pivot.columns.droplevel()
    pivot['sum']=pivot.sum(axis=1)
    pivot=pivot.sort_values(by='sum',ascending=False)
    pivot=pivot.drop('sum',axis=1)
    pivot.to_csv(f'{tab_path}/{sample}_matrix.tsv',sep='\t')
    if len(pivot)>10:
        pivot=pivot[0:10]
    g=sns.heatmap(pivot,square=square,cbar_kws={'fraction' : 0.01},cmap='OrRd',linewidth=1,annot=annot,fmt='.1f')
    return g


with PdfPages(snakemake.output.heatmap_counts) as pdf:
    for sample in unique_sample_list:
        plt.figure(figsize=(11.7,8.27))
        plt.subplots_adjust(left=0.3, bottom=0.3)
        plt.title(sample)
        fig=draw_heatmap_counts(all_tools,sample,rank,square,annot,path)
        pdf.savefig()
        plt.close()



def get_presence_table(table, rank, tool, sample):
    table=table[table['sample']==sample]
    all_found=list(set(table[f'{rank}']))
    tool_found=list(set(table[table['tool']==f'{tool}'][f'{rank}']))
    dic={}
    for hit in all_found:
        dic[hit]={}
        if hit in tool_found:
            dic[hit][f'{tool}'] = True
        else:
            dic[hit][f'{tool}'] = False
    df=pd.DataFrame.from_dict(dic, orient='index')
    df.index.name=f'{rank}'
    return df

def draw_presence_heatmap_tools(table,rank,sample,fraction,squared):
    tool_list=list(set(table['tool']))
    tb_list=[]
    for tool in tool_list:
        tb_list.append(get_presence_table(all_tools,rank,tool,sample))
    matrix=pd.concat(tb_list,axis=1,sort=False)
    matrix['sum']=matrix.sum(axis=1)
    matrix=matrix.sort_values(by='sum',ascending=False)
    max_hits_per_tool=max(matrix['sum'])
    frac=fraction
    min_hits_to_consider=int(np.floor(max_hits_per_tool*frac))
    subset=matrix[matrix['sum'].isin(range(min_hits_to_consider, max_hits_per_tool+1))]#Consider only the genomes identified by 4 out of 7 tools
    subset=subset.drop('sum',axis=1)
    categories=sorted(np.unique(matrix))
    if len(categories)==1 and categories[0]==True:
        mycolors = [(0.3, 0.6, 0.1, 1.0),(0.0, 0.0, 0.0, 0.9)]
    else:
        mycolors = [(0.0, 0.0, 0.0, 0.9), (0.3, 0.6, 0.1, 1.0)]#RGB color palette
    cmap = LinearSegmentedColormap.from_list('Custom', mycolors, len(mycolors))
    if len(subset)>=10:#limit matrix length to 10
        subset=subset[0:10]
    g = sns.heatmap(subset, square=squared, cbar_kws={'fraction' : 0.01},cmap=cmap,linecolor='white',linewidths='0.1')
    colorbar = g.collections[0].colorbar
    colorbar.set_ticks([0,1.0])
    colorbar.set_ticklabels(['absent', 'present'])
    return g

#Draw heatmap presence between tools for the same samples
with PdfPages(snakemake.output.heatmap_presence_tools) as pdf:
    for sample in unique_sample_list:
        plt.figure(figsize=(11.7,8.27))
        plt.subplots_adjust(left=0.3, bottom=0.3)
        plt.title(sample)
        fig=draw_presence_heatmap_tools(all_tools,rank,sample,fraction,square)
        pdf.savefig()
        plt.close()


def draw_heatmap_presence_samples(table, rank, tool):
    t = table[table['tool'] == tool]
    ranks_list = list(t[f'{rank}'])
    hits = sorted(set(ranks_list), key=ranks_list.index)  # Keep ranks sorted by original order (read counts)
    samples = list(set(t['sample']))
    dic = {}
    for hit in hits:
        dic[hit]={}
        for s in samples:
            hits_found_sample=list(set(t[t['sample'] == f'{s}'][f'{rank}']))
            if hit in hits_found_sample:
                dic[hit][s] = True
            else:
                dic[hit][s] = False
    df = pd.DataFrame.from_dict(dic, orient='index')
    df.index.name = f'{rank}'
    categories = sorted(np.unique(df))
    if len(categories) == 1 and categories[0] == True:
        mycolors = [(0.3, 0.6, 0.1, 1.0), (0.0, 0.0, 0.0, 0.9)]
    else:
        mycolors = [(0.0, 0.0, 0.0, 0.9), (0.3, 0.6, 0.1, 1.0)]  # RGB color palette
    cmap = LinearSegmentedColormap.from_list('Custom', mycolors, len(mycolors))
    hm = sns.heatmap(df, square=False, cbar_kws={'fraction': 0.01}, cmap=cmap, linecolor='white', linewidths='0.1')
    colorbar = hm.collections[0].colorbar
    colorbar.set_ticks([0, 1.0])
    colorbar.set_ticklabels(['absent', 'present'])
    return hm


tools=list(set(all_tools['tool']))

#Draw heatmap presence between samples for the same tool
with PdfPages(snakemake.output.heatmap_presence_samples) as pdf:
    for tool in tools:
        plt.figure(figsize=(11.7, 8.27))
        plt.subplots_adjust(left=0.3, bottom=0.3)
        plt.title(tool)
        draw_heatmap_presence_samples(all_tools, rank, tool)
        pdf.savefig()
        plt.close()
