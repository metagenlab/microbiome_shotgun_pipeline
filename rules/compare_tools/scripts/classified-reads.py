import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



superkingdom=snakemake.params.superkingdom
def groupby_samples(table,tool_name,rank):
    if tool_name=='ezvir':
        sample_tb=table.groupby(['sample',f'{rank}'],as_index=False)['read_counts'].sum()
        sample_tb['read_counts']=sample_tb['read_counts']/2
    else:
        sample_tb = table.groupby(['sample', f'{rank}'], as_index=False)['read_counts'].sum()
    sample_tb['tool'] = [tool_name] * len(sample_tb)
    sample_tb['rank']=[rank]*len(sample_tb)
    return sample_tb

rank='superkingdom'
tables=[]
files=snakemake.input.all_tools
for file in files:
    toolname = file.split('/')[1]
    table=pd.read_csv(file,sep='\t')
    table=table[table.superkingdom==superkingdom]
    gbs=[]
    tb=groupby_samples(table,toolname,rank)
    gbs.append(tb)
    gbs=pd.concat(gbs,sort=False)
    tables.append(gbs)
tables=pd.concat(tables,sort=False)
tables=tables.groupby(['tool','sample','rank'],as_index=False)['read_counts'].sum()
tables.to_csv(snakemake.output.table,sep='\t',index=False)

gold_standard=pd.read_csv(snakemake.input.gold_standard,sep='\t')
gold_standard=gold_standard[gold_standard.superkingdom==superkingdom]
gold_standard=gold_standard.groupby(['sample',f'{rank}'],as_index=False).sum()
theoretical_max=max(gold_standard['read_counts'])
plt.figure(figsize=(11.7,8.27))
plt.subplots_adjust(left=0.4,right=0.81)
bp=sns.barplot(data=tables,y='read_counts',x='tool')
ax1=bp.axes
ax1.axhline(theoretical_max,ls='-',c='red')
ax1.text(1,theoretical_max+1000,'theoretical read counts',c='red')
bp.get_figure().savefig(snakemake.output.barplot)