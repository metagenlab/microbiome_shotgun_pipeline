import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

rank=snakemake.params["rank"]
tb=pd.read_csv(snakemake.input[0],sep='\t')
tb=tb.replace({"TP":0,"FP":2,"FN":1})
subset=tb[[f'{rank}','tool','present']]

matrix=pd.pivot(data=subset,index=f'{rank}',columns='tool',values='present')
matrix=matrix.replace(np.nan,3)
col_list=list(matrix.columns)
matrix=matrix.sort_values(by=col_list,axis=0,ascending=True)
if len(matrix)>50:
    matrix=matrix[0:50]

if 1 not in np.unique(matrix):#if no false negatives then no black colour
    myColors = [(0.3, 0.6, 0.1, 1.0),(1.0, 0.64, 0, 1.0),(0, 0, 0, 0)] #(Green, Orange, White)
if 2 not in np.unique(matrix):#if no false positives then no orange
    myColors=[(0.3, 0.6, 0.1, 1.0),(0.0, 0.0, 0.0, 0.9),(0, 0, 0, 0)]
if 3 not in np.unique(matrix):
    myColors = [(0.3, 0.6, 0.1, 1.0),(0.0, 0.0, 0.0, 0.9), (1.0, 0.64, 0, 1.0)]
else:
    myColors = [(0.3, 0.6, 0.1, 1.0),(0.0, 0.0, 0.0, 0.9), (1.0, 0.64, 0, 1.0),
                (0, 0, 0, 0)]#Green Black orange white

cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

plt.figure(figsize=(11.7,8.27))
sns.set(font_scale=0.6)
plt.subplots_adjust(left=0.8)
ax=sns.heatmap(matrix,linecolor='black',linewidths='0.1',cmap=cmap,cbar_kws={'fraction' : 0.04},square=True)
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0,1.0,2.0,3.0])
colorbar.set_ticklabels(["TP","FN","FP","NA"])
ax.get_figure().savefig(snakemake.output[0])