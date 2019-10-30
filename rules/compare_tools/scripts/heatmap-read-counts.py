
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

rank=snakemake.params["rank"]
column=snakemake.params["value"]#select read counts or read percentages

tb=pd.read_csv(snakemake.input[0],sep='\t')
col=[]
for colname in tb.columns:
    if column in colname:
        col.append(colname)

subset=tb.groupby(rank).sum()
subset=subset[col]

subset=subset.sort_values(by=col,ascending=False)

plt.figure(figsize=(20,10))
sns.set(font_scale=0.5)
hm=sns.heatmap(subset,cmap="YlGnBu",fmt=".1f")

hm.get_figure().savefig(snakemake.output[0])
