import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

rank=snakemake.params["rank"]
tb=pd.read_csv(snakemake.input[0],sep='\t')
tb=tb.replace({"TP":0,"FP":2,"FN":1})
subset=tb[['tool','sample',f'{rank}','present']]

def draw_hm(subset,sample,rank):
    subset=subset[subset['sample']==sample]
    matrix=pd.pivot(data=subset,index=f'{rank}',columns='tool',values='present')
    matrix=matrix.replace(np.nan,3)
    col_list=list(matrix.columns)
    matrix=matrix.sort_values(by=col_list,axis=0,ascending=True)
    if len(matrix)>30:
        matrix=matrix[0:30]

    if 1 not in np.unique(matrix):#if no false negatives then no black colour
        myColors = [(0.3, 0.6, 0.1, 1.0),(1.0, 0.64, 0, 1.0),(0, 0, 0, 0)] #(Green, Orange, White)
    if 2 not in np.unique(matrix) and 3 not in np.unique(matrix):#if no false positives then no orange
        myColors=[(0.3, 0.6, 0.1, 1.0),(0.0, 0.0, 0.0, 0.9)]
    else:
        myColors = [(0.3, 0.6, 0.1, 1.0),(0.0, 0.0, 0.0, 0.9), (1.0, 0.64, 0, 1.0),
                    (0, 0, 0, 0)]#Green Black orange white

    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    plt.figure(figsize=(11.7,8.27))
    sns.set(font_scale=0.6)
    plt.subplots_adjust(left=0.8)
    plt.title(f'sample={sample}')
    ax=sns.heatmap(matrix,linecolor='black',linewidths='0.1',cmap=cmap,cbar_kws={'fraction' : 0.04},square=True)
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0,1.0,2.0,3.0])
    colorbar.set_ticklabels(["TP","FN","FP","NA"])
    return ax

samples=list(set(subset['sample']))
with PdfPages(snakemake.output[0]) as pdf:
    for sample in samples:
        fig=draw_hm(subset,sample,rank)
        pdf.savefig()
        plt.close()
