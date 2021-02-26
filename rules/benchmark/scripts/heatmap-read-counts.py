import pandas as pd
import seaborn as sns
from fpdf import FPDF
import matplotlib.pyplot as plt


rank=snakemake.params["rank"]
values=snakemake.params["value"]#select read counts or read percentages
superkingdom=snakemake.params["superkingdom"]

def get_heatmap_table(file,rank,value,superkingdom):
    tb=pd.read_csv(file,sep='\t')
    subset = tb[tb.superkingdom.str.lower() == superkingdom]
    subset=subset.groupby(f'{rank}',as_index=False).sum()
    subset.set_index(f'{rank}',inplace=True)
    col=[]
    for colname in subset.columns:
        if value in colname:
            col.append(colname)
    gb=subset[col]
    gb=gb.sort_values(by=col,axis=0,ascending=False)
    return gb



tb=get_heatmap_table(snakemake.input[0],rank,values,superkingdom)



if len(tb)>0:
    if len(tb)>10:
        tb=tb[0:10]
        plt.figure(figsize=(11.7,8.27))
        plt.subplots_adjust(left=0.4, bottom=0.4)
        sns.set(font_scale=1.0)
        hm=sns.heatmap(tb,cmap="YlGnBu")
        plt.xticks(rotation=90)
        hm.get_figure().savefig(snakemake.output[0])
else:
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt=f"No {superkingdom} found", ln=1, align="C")
    pdf.output(snakemake.output[0])



