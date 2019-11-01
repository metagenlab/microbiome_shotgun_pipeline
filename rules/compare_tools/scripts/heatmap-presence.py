import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
rank=snakemake.params["rank"]
tb=pd.read_csv(snakemake.input[0],sep='\t')
tb=tb.replace(True,1)
subset=tb[[f'{rank}','tool','present']]
matrix=pd.pivot_table(data=subset,index=subset['species'],columns=subset['tool'])
matrix.columns=matrix.columns.droplevel()
plt.figure(figsize=(11.7,8.27))
sns.set(font_scale=1.0)
myColors = ((0.0, 0.0, 0.0, 0.9), (0.3, 0.6, 0.1, 1.0))#RGB color palette
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
plt.subplots_adjust(left=0.2,bottom=0.2)
ax=sns.heatmap(matrix,linecolor='black',linewidths='0.1',cmap=cmap,square=True)
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0,1.0])
colorbar.set_ticklabels(['absent', 'present'])
ax.get_figure().savefig(snakemake.output[0])