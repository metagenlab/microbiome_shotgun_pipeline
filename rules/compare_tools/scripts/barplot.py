import pandas as pd
import seaborn as sns
from numpy import mean
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
tool_name=snakemake.wildcards.tool
rank=snakemake.params.rank

superkingdoms=['Bacteria','Viruses']



gold_standard = pd.read_csv(snakemake.input.gold_standard,sep='\t',index_col='taxid')
gold_standard=gold_standard.drop(9606)#drop human
gold_standard=gold_standard.reset_index()
gold_standard = gold_standard.groupby([f'{rank}','superkingdom'],as_index=False)['read_counts'].sum()
gold_standard['type'] = ['gold_standard'] * len(gold_standard)# add a column to specify if the read counts are from a tool or the gold standard

tool_out = pd.read_csv(snakemake.input.tool_out,sep='\t')



tool_out=tool_out.groupby([f'{rank}','superkingdom'],as_index=False)['read_counts'].sum()
tool_out['type']=[tool_name]*len(tool_out)
cat = pd.concat([gold_standard, tool_out], axis=0, sort=False)

cat.set_index(f'{rank}',inplace=True)
true_genomes_list=list(cat.loc[cat['type']=='gold_standard'].index)
cat=cat.reset_index()
subset=cat.loc[cat[f'{rank}'].isin(true_genomes_list)]
subset=subset.sort_values(by='read_counts',ascending=False)




def draw_barplot(table,superkingdom,tool_name,estimator):
    plt.figure(figsize=(11.7,8.27))
    plt.subplots_adjust(left=0.4,right=0.81)
    plt.title(superkingdom)
    sns.set(font_scale=0.9)
    table=table[table.superkingdom==superkingdom]
    palette = {'gold_standard': 'C0', f'{tool_name}': 'C1'}
    bp=sns.barplot(data=table,y=f'{rank}',x='read_counts',hue='type',palette=palette,estimator=estimator)
    bp.legend(loc=2, bbox_to_anchor=(1.05, 1),borderaxespad=0.1)
    return bp

with PdfPages(snakemake.output[0]) as pdf:
    for superkingdom in superkingdoms:
        fig=draw_barplot(subset,superkingdom,tool_name,mean)
        pdf.savefig()
        plt.close()


