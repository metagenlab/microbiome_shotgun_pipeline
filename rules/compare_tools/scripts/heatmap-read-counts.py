import pandas as pd
import seaborn as sns
from fpdf import FPDF
import matplotlib.pyplot as plt


rank=snakemake.params["rank"]
values=snakemake.params["value"]#select read counts or read percentages
superkingdom=snakemake.params["superkingdom"]

def get_heatmap_table(file,rank,value,superkingdom):
    tb=pd.read_csv(file,sep='\t')
    col=[]
    for colname in tb.columns:
        if value in colname:
            col.append(colname)

    subset = tb[tb.superkingdom == superkingdom]

    gb=subset.groupby(rank).sum()

    gb=gb[col]
    gb=gb.sort_values(by=col,ascending=False)
    return gb



tb=get_heatmap_table(snakemake.input[0],rank,values,superkingdom)

if len(tb)>0:
    plt.figure(figsize=(21, 10))
    sns.set(font_scale=0.5)
    hm=sns.heatmap(tb,cmap="YlGnBu")
    plt.xticks(rotation=45)
    hm.get_figure().savefig(snakemake.output[0])
else:
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt=f"No {superkingdom} found", ln=1, align="C")
    pdf.output(snakemake.output[0])



