import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
tables=snakemake.input.all_tools
rank=snakemake.params.rank
superkingdom=snakemake.params.superkingdom
gold_standard_file=snakemake.input.gold_standard

gold_standard=pd.read_csv(gold_standard_file,sep='\t')
true_superkingdom=gold_standard[gold_standard['superkingdom']==superkingdom]
true_superkingdom=true_superkingdom.groupby([f'{rank}','sample'],as_index=False).sum()


def get_l1_distance(true_tb,tool_tb,tool_name,rank):
    sample_list=sorted(list(set(tool_tb['sample'])))
    gs_sample_list=sorted(list(set(true_tb['sample'])))

    tables=[]
    for sample in sample_list:
        sample_tb=tool_tb[tool_tb['sample']==sample]
        if len(gs_sample_list)==1:
            true_tb=true_tb
        else:
            true_tb=true_tb[true_tb['sample']==sample]
        sample_tb=sample_tb.groupby(f'{rank}',as_index=False).sum()
        if tool_name=='ezvir':
            sample_tb['read_counts']=sample_tb['read_counts']//2
        ljoin=pd.merge(true_tb,sample_tb,how='left',on=f'{rank}',suffixes=('_true','_found'))
        #ljoin.replace(np.nan,0,inplace=True)
        ljoin['l1_distance']=abs(ljoin['read_counts_true']-ljoin['read_counts_found'])
        ljoin['tool']=[tool_name]*len(ljoin)
        ljoin['sample']=[sample]*len(ljoin)
        subset=ljoin[[f'{rank}','read_counts_true','read_counts_found','l1_distance','sample','tool']]
        tables.append(subset)
    ljoin_all=pd.concat(tables,sort=False)
    return ljoin_all

tool_type=snakemake.params.type
all_tools_l1=[]
for file in tables:
    tool = os.path.split(file)[1].split('.tsv')[0]
    type = tool_type[tool]
    tb=pd.read_csv(file,sep='\t')
    tool_tb=tb[tb.superkingdom==superkingdom]
    tool_l1 = get_l1_distance(true_superkingdom, tool_tb, tool, rank)
    tool_l1['type'] = [type] * len(tool_l1)
    all_tools_l1.append(tool_l1)

all_tools=pd.concat(all_tools_l1,sort=False)
all_tools=all_tools.dropna()
all_tools.reset_index(drop=True,inplace=True)
all_tools.to_csv(snakemake.output.l1_dist_table,sep='\t',index=False)
plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(right=0.81,bottom = 0.3)
sns.set_context('paper',font_scale=1.5)
vp=sns.violinplot(data=all_tools,x='tool',y='l1_distance',hue='type',dodge=False,cut=0)
vp.set_xticklabels(vp.get_xticklabels(),rotation=90)
vp.legend(loc=2, bbox_to_anchor=(1.05, 1),borderaxespad=0.1)
vp.get_figure().savefig(snakemake.output.l1_distance)
