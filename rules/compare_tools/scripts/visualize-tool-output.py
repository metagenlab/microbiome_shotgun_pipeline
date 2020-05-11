import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from numpy import mean
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

def concat_sample_counts(file,target_val,sample_name, rank, tool_name, superkingdom, threshold):
    tb_list = []
    col_sel = []
    tool_tb = pd.read_csv(file, sep='\t')
    tool_tb = tool_tb[tool_tb.superkingdom == superkingdom]
    tool_tb = tool_tb.groupby(f'{rank}').sum()

    for col in tool_tb:
        if target_val in col:
            col_sel.append(col)
    subset = tool_tb[col_sel]

    for col in subset.columns:
        tb = subset[[col]]
        tb = tb[tb > threshold].dropna()
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
bp_threshold=snakemake.params.threshold
superkingdom=snakemake.params.superkingdom
threshold=0
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
if value=="counts":
    bp_table=all_tools[all_tools.read_counts>=bp_threshold].dropna()
else:
    bp_table = all_tools[all_tools.read_percent >= bp_threshold].dropna()


unique_sample_list=list(set(bp_table['sample']))
unique_names=list(set(bp_table[f'{rank}']))

fontsize=10
square=False
annot=True
fraction=0.6#For the presence heatmap: only consider hits identified by 60% of the tools
if len(unique_names)>=150:
    fontsize=2
    square=False
    fraction=0.9#High number of hits = heatmap is not visible
if len(unique_names)<=150 and len(unique_names)>100:
    fontsize=3
    square=False
elif len(unique_names)<100:
    fontsize=7
    square=True




sns.set(font_scale=1.0)
plt.figure(figsize=(len(bp_table)/2,len(bp_table)))
plt.subplots_adjust(left=0.2,top=0.8)
bp=sns.catplot(data=bp_table,x=f'read_{value}',y=f'{rank}',col='sample',col_wrap=3,kind='bar',estimator=mean,errwidth=0.3)
bp.set(xscale='log')

for ax in bp.axes.flat:
    for ylab in ax.get_yticklabels():
        ylab.set_size(fontsize)
    for xlab in ax.get_xticklabels():
        xlab.set_size(10)

bp.savefig(snakemake.output.barplot_samples)


sns.set(font_scale=1.0)
plt.figure(figsize=(len(bp_table)/2,len(bp_table)))
plt.subplots_adjust(left=0.2,top=0.8)
bpt=sns.catplot(data=bp_table,x=f'read_{value}',y=f'{rank}',hue='sample',col='tool',col_wrap=3,kind='bar',estimator=mean,errwidth=0.3)
bpt.set(xscale='log')

for ax in bpt.axes.flat:
    for ylab in ax.get_yticklabels():
        ylab.set_size(fontsize)
    for xlab in ax.get_xticklabels():
        xlab.set_size(10)

bpt.savefig(snakemake.output.barplot_tools)




def draw_heatmap_counts(data,sample,rank,square,annot):
    subset=data[data['sample']==sample]
    pivot=subset.pivot_table(subset,index=f'{rank}',columns='tool')
    pivot=pivot.replace(np.nan,0)
    pivot.columns=pivot.columns.droplevel()
    pivot['sum']=pivot.sum(axis=1)
    pivot=pivot.sort_values(by='sum',ascending=False)
    pivot=pivot.drop('sum',axis=1)
    #pivot=pivot.replace(0,1)
    #logp=np.log10(pivot)
    g=sns.heatmap(pivot,square=square,cbar_kws={'fraction' : 0.01},cmap='OrRd',linewidth=1,annot=annot,fmt='.1f')
    return g

fonstscale=1.0
if len(bp_table)>100:
    fonstscale=5.0
if len(bp_table)>=50 and len(bp_table)<=100:
    fonstscale=2.0

with PdfPages(snakemake.output.heatmap_counts) as pdf:
    for sample in unique_sample_list:
        plt.figure(figsize=(len(bp_table)/2,len(bp_table)))
        sns.set(font_scale=fonstscale)
        plt.subplots_adjust(left=0.3, bottom=0.3)
        plt.title(sample)
        fig=draw_heatmap_counts(bp_table,sample,rank,square,annot)
        pdf.savefig()
        plt.close()



def get_presence_table(table,rank,tool,sample):
    table=table[table['sample']==sample]
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

def draw_presence_heatmap(table,rank,sample,fraction,squared):
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
    subset=matrix[matrix['sum'].isin(range(min_hits_to_consider,max_hits_per_tool+1))]#Consider only the genomes identified by 4 out of 7 tools
    subset=subset.drop('sum',axis=1)
    categories=sorted(np.unique(matrix))
    if len(categories)==1 and categories[0]==True:
        myColors = [(0.3, 0.6, 0.1, 1.0),(0.0, 0.0, 0.0, 0.9)]
    else:
        myColors = [(0.0, 0.0, 0.0, 0.9), (0.3, 0.6, 0.1, 1.0)]#RGB color palette
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
    g = sns.heatmap(subset, square=squared, cbar_kws={'fraction' : 0.01},cmap=cmap,linecolor='white',linewidths='0.1')
    colorbar = g.collections[0].colorbar
    colorbar.set_ticks([0,1.0])
    colorbar.set_ticklabels(['absent', 'present'])
    return g


with PdfPages(snakemake.output.heatmap_presence) as pdf:
    for sample in unique_sample_list:
        plt.figure(figsize=(len(bp_table),len(bp_table)/2))
        sns.set(font_scale=fonstscale)
        plt.subplots_adjust(left=0.3, bottom=0.3)
        plt.title(sample)
        fig=draw_presence_heatmap(all_tools,rank,sample,fraction,square)
        pdf.savefig()
        plt.close()



