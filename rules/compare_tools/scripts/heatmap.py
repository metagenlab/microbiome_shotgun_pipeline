
import pandas as pd
import seaborn as sns

rank=snakemake.params["rank"]


tb=pd.read_csv(snakemake.input[0],sep='\t')
col=[]
for colname in tb.columns:
    if 'counts' in colname:
        col.append(colname)

subset=tb.groupby(rank).sum()
subset=subset[col]

subset=subset.sort_values(by=col,ascending=False)
plot=sns.heatmap(subset,cmap="YlGnBu")
fig=plot.get_figure()
fig.savefig(snakemake.output[0])