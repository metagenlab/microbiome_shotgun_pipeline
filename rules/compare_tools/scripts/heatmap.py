
import pandas as pd
import seaborn as sns

rank=snakemake.params["rank"]


tb=pd.read_csv(snakemake.input[0],sep='\t')
subset=tb.groupby([rank]).sum()
subset=subset.sort_values(axis=1,ascending=False)
heatmap=sns.heatmap(subset,cmap="YlGnBu")
heatmap.savefig(smakemake.output[0])

