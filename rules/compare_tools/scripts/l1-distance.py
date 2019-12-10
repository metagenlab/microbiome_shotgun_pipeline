import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

tables=snakemake.input.all_tools
rank=snakemake.params.rank
superkingdom=snakemake.params.superkingdom
gold_standard_file=snakemake.input.gold_standard

gold_standard=pd.read_csv(gold_standard_file,sep='\t')
true_superkingdom=gold_standard[gold_standard['superkingdom']==superkingdom]



def get_l1_distance(true_tb,tool_tb,tool_name,rank):
    sample_list=sorted(list(set(tool_tb['sample'])))
    tables=[]
    for sample in sample_list:
        sample_tb=tool_tb[tool_tb['sample']==sample]
        sample_tb=sample_tb.groupby(f'{rank}',as_index=False).sum()
        ljoin=pd.merge(true_tb,sample_tb,how='left',on=f'{rank}',suffixes=('_true','_found'))
        #ljoin.replace(np.nan,0,inplace=True)
        ljoin['l1_distance']=ljoin['read_counts_true']-ljoin['read_counts_found']
        ljoin['tool']=[tool_name]*len(ljoin)
        ljoin['sample']=[sample]*len(ljoin)
        subset=ljoin[[f'{rank}','read_counts_true','read_counts_found','l1_distance','sample','tool']]
        tables.append(subset)
    ljoin_all=pd.concat(tables,sort=False)
    return ljoin_all

tool_type={'surpi':'mapping','pathseq':'mapping','ezvir':'mapping','ganon':'k-mer','kraken2':'k-mer','kaiju':'k-mer','bracken':'k-mer'}
all_tools_l1=[]
for file in tables:
    tool = file.split('/')[1]
    tb=pd.read_csv(file,sep='\t')
    tool_tb=tb[tb.superkingdom==superkingdom]
    tool_l1 = get_l1_distance(true_superkingdom, tool_tb, tool, rank)
    tool_l1['type'] = [tool_type[tool]] * len(tool_l1)
    all_tools_l1.append(tool_l1)

all_tools=pd.concat(all_tools_l1,sort=False)
all_tools=all_tools.dropna()
all_tools.reset_index(drop=True,inplace=True)
all_tools.to_csv(snakemake.output.l1_dist_table,sep='\t',index=False)
plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(right=0.81)
vp=sns.violinplot(data=all_tools,x='tool',y='l1_distance',hue='type',dodge=False)
vp.legend(loc=2, bbox_to_anchor=(1.05, 1),borderaxespad=0.1)
vp.get_figure().savefig(snakemake.output.l1_distance)
