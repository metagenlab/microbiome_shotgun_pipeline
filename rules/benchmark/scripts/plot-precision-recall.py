import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

all_tb = []
files = snakemake.input
for file in files:
    tb = pd.read_csv(file,sep='\t')
    all_tb.append(tb)

concat = pd.concat(all_tb)
sorted_aupr = concat.sort_values(by=['AUPR','tool','sample'],ascending=False)
sorted_aupr.to_csv(snakemake.output.curve_table,sep='\t',index=False)

plt.figure(figsize=(11.7,8.27))
sns.set_context('paper',font_scale=1.2)
prc=sns.FacetGrid(data=sorted_aupr,col='tool',col_wrap=3,hue='AUPR')
prc.map(plt.step,'recall','precision',where='pre',label='AUPR')
prc.map(plt.fill_between,'recall','precision',alpha=0.4,step='pre')
prc.add_legend()
prc.savefig(snakemake.output.pr_curve)

plt.figure(figsize=(11.7,8.27))
sns.set(context='paper',style='darkgrid',font_scale=1.2)
f1_plot=sns.FacetGrid(data=sorted_aupr,col='tool',col_wrap=3)
f1_plot.map(plt.plot,'threshold','precision',label='precision',color='blue',ls=':',alpha=0.6)
f1_plot.map(plt.plot,'threshold','recall',label='recall',color='red',ls=':',alpha=0.6)
f1_plot.map(plt.plot,'threshold','F1_score',label='F1',color='green')
f1_plot.set_axis_labels('threshold','score')
for ax in f1_plot.axes:
    tool_name=ax.title.get_text().split('tool = ')[1]
    max_f1=max(sorted_aupr[sorted_aupr.tool == tool_name].F1_score)
    index_max=sorted_aupr[sorted_aupr.tool == tool_name]
    min_cutoff=list(index_max[index_max.F1_score==max_f1].threshold)[0]
    ax.axvline(min_cutoff,ls='-',color='black',label='optimal-threshold',alpha=0.4)
    ax.text(min_cutoff+300,0,f'{min_cutoff}',color='black',alpha=0.4)
f1_plot.add_legend()
f1_plot.savefig(snakemake.output.f1_curve)

indexed_subset=sorted_aupr.set_index('threshold')
opt_read_thresh={}

tool_list=list(set(sorted_aupr['tool']))
for tool in tool_list:
    tool_tab=indexed_subset[indexed_subset.tool==tool]
    max_score=max(tool_tab.F1_score)
    index_list=list(tool_tab[tool_tab.F1_score==max_score].index)
    optimal_threshold=index_list[0]
    opt_read_thresh[tool]={'read_threshold':optimal_threshold,'max_F1_score':max_score}

opt_thresholds_tb=pd.DataFrame.from_dict(opt_read_thresh,orient='index')
opt_thresholds_tb.to_csv(snakemake.output.optimal_params,sep='\t')