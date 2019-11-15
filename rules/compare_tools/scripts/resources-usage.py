import pandas as pd
import math
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import mean
resources=snakemake.input

all_tools=[]
tool_type={'surpi':'mapping','pathseq':'mapping','ezvir':'mapping','ganon':'k-mer','kraken2':'k-mer','kaiju':'k-mer','bracken':'k-mer'}
for file in resources:
    tool=file.split('/')[3].split('.')[0]
    name=file.split('/')[2]
    tb = pd.read_csv(file,sep='\t')
    tb['tool']=[tool]*len(tb)
    tb['type']=tool_type[tool]
    tb['sample']=[name]*len(tb)
    all_tools.append(tb)

all_tools=all_tools.replace('-',0)


all_tools['max_mem_gb']=all_tools['max_rss']/math.exp(9)
all_tools['minutes']=all_tools['s']/60

all_tools.to_csv(snakemake.output.tab)

plt.figure(figsize=(11.7,8.27))
bpm=sns.catplot(data=all_tools,y='tool',x='max_mem_gb',col='type',kind='bar',estimator=mean)
bpm.set(xlabel='GB')
for ax in bpm.axes.flat:
    ax.set_ylabel(ax.get_ylabel(), rotation=0)
bpm.get_figure().savefig(snakemake.output.memo)

plt.figure(figsize=(11.7,8.27))
bpt=sns.catplot(data=all_tools,y='tool',x='minutes',col='type',kind='bar',estimator=mean)
bpt.set(xlabel='minutes')
for ax in bpt.axes.flat:
    ax.set_ylabel(ax.get_ylabel(), rotation=0)
bpt.get_figure().savefig(snakemake.output.time)