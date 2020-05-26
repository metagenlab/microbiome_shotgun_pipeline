import pandas as pd
from Bio import Entrez
import matplotlib.pyplot as plt

Entrez.email=snakemake.params.NCBI_email
Entrez.api_key=snakemake.params.NCBI_key

samplename=snakemake.wildcards.sample
fullname=snakemake.wildcards.filename

temp=fullname.split('_')
acc='_'.join(temp[:2])#split on the second underscore to get accession id
ID=Entrez.read(Entrez.esearch(db='assembly', term=f'{acc}'))['IdList'][0]
summary=Entrez.read(Entrez.esummary(db='assembly', id=ID))['DocumentSummarySet']['DocumentSummary'][0]
name=summary['Organism']
species=summary['SpeciesName']
taxid=summary['Taxid']
accession=summary['AssemblyAccession']

coverage=pd.read_csv(snakemake.input.table,sep='\t',names=['header','pos','cov'])
f=open(snakemake.input.txt,'r')
info=f.read()
reads_mapped=info.split('\n')[4].split(' ')[0]

covered=len(coverage[coverage['cov']!=0])#Coverage = positions that are mapped by at least one read
covpercent=round((covered/len(coverage))*100,2)#Coverage percent = Coverage/genome length *100 (rounded to 2 decimals)

plt.figure(figsize=(11.7,8.27))
plt.title(f'{samplename} reads mapped to {name}, assembly accession: {accession}\n\n {covpercent}% covered\n\n total reads mapped:{reads_mapped}',fontsize=12,ha='center')
plt.plot('pos','cov',data=coverage)
plt.xlabel('position',fontsize=11)
plt.ylabel('depth',fontsize=11)
plt.savefig(snakemake.output[0])
