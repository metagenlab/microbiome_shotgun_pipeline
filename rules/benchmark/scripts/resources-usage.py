import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import mean
import numpy as np
resources=snakemake.input

tb_list=[]
tool_type=snakemake.params.type
for file in resources:
    tool=file.split('/')[3].split('.')[0]
    name=file.split('/')[2]
    tb = pd.read_csv(file,sep='\t')
    tb['tool']=[tool]*len(tb)
    tb['type']=tool_type[tool]
    tb['sample']=[name]*len(tb)
    tb_list.append(tb)

all_tools=pd.concat(tb_list)
all_tools=all_tools.replace('-',0)


all_tools['max_mem_gb']=all_tools['max_rss']/1024
all_tools['minutes']=all_tools['s']/60
all_tools['n_cpu_usage']=all_tools['mean_load']/100

all_tools.to_csv(snakemake.output.tab,sep='\t')

plt.figure(figsize=(11.7,8.27))
bpm=sns.catplot(data=all_tools,x='tool',y='max_mem_gb',col='type',kind='bar',estimator=mean,log=True)
bpm.set(ylabel='GB')
for ax in bpm.axes.flat:
    for label in ax.get_xticklabels():
        label.set_rotation(90)
bpm.savefig(snakemake.output.memo)

plt.figure(figsize=(11.7,8.27))
bpt=sns.catplot(data=all_tools,x='tool',y='minutes',col='type',kind='bar',estimator=mean)
bpt.set(ylabel='minutes')
for ax in bpt.axes.flat:
    for label in ax.get_xticklabels():
        label.set_rotation(90)
bpt.savefig(snakemake.output.time)

fig = plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(bottom = 0.3)
plt.rcParams.update({'font.size': 14})
ax = fig.add_subplot(111)
ax2 = ax.twinx()
width = 0.4
gb=all_tools.groupby('tool',as_index=False).mean()
N = gb.shape[0]
dim = 2
lim = (dim - 1) * 0.4
offsets = np.linspace(-lim, lim, dim)
ps=[]
ps.append(ax.bar(data=gb,x=np.arange(N)*dim + offsets[0],height='max_mem_gb',color='salmon',label='memory (GB)'))
ps.append(ax2.bar(data=gb,x=np.arange(N)*dim + offsets[1],height='minutes',color='pink',label='minutes'))
ax.set_xticks(np.arange(N) * dim)
ax.set_xticklabels([tool for tool in gb['tool']],rotation=90)
ax.legend(ps,['memory','time'],loc='upper left')
ax.set_ylabel('Gigabytes')
ax2.set_ylabel('Minutes')
fig.savefig(snakemake.output.mem_and_time)