import pandas as pd
from Bio import Entrez
import matplotlib.pyplot as plt

Entrez.email=snakemake.params.NCBI_email
Entrez.api_key=snakemake.params.NCBI_key

samplename=snakemake.wildcards.sample
fullname=snakemake.wildcards.filename


temp=fullname.split('_')
acc='_'.join(temp[:2])#split on the second underscore to get accession id
try:
    ID=Entrez.read(Entrez.esearch(db='assembly', term=f'{acc}'),validate=False)['IdList'][0]
    summary=Entrez.read(Entrez.esummary(db='assembly', id=ID),validate=False)['DocumentSummarySet']['DocumentSummary'][0]
    name=summary['Organism']
    species=summary['SpeciesName']
    taxid=summary['Taxid']
    accession=summary['AssemblyAccession']
except IndexError:
    name=fullname
    species='unknown'
    taxid='unknown'
    accession='unknown'



coverage=pd.read_csv(snakemake.input.table,sep='\t',names=['header','pos','cov'])
f=open(snakemake.input.txt,'r')
info=f.read()
reads_mapped=info.split('\n')[4].split(' ')[0]

covered=len(coverage[coverage['cov']!=0])#Coverage = positions that are mapped by at least one read
covpercent=round((covered/len(coverage))*100,2)#Coverage percent = Coverage/genome length *100 (rounded to 2 decimals)
med_reads=round(coverage['cov'].median(),2)#Calculate median read coverage
rollm=coverage.rolling(len(coverage)//2).mean()#calculate moving average with a sliding window of genome length divided by 100

plt.style.use('ggplot')
plt.figure(figsize=(11.7,8.27))
plt.title(f'{samplename} reads mapped to {name}, assembly accession: {accession}\n\n{covpercent}% covered',fontsize=14,ha='center')
plt.gcf().text(0.92, 0.5,f'Total reads mapped: {reads_mapped}\n\nMedian reads mapped: {med_reads}',fontsize=12,bbox={'facecolor':'white','alpha':0.7,'edgecolor': 'black','pad':0.5,'boxstyle':'round'},style='italic')
plt.plot('pos','cov',data=coverage,c='darkblue',alpha=0.2,label='read number')
plt.plot('pos','cov',data=rollm,c='blue',label='rolling mean')
plt.legend(loc=2, bbox_to_anchor=(1.05, 1),borderaxespad=0.1)
plt.ylabel('depth')
plt.xlabel('position')
plt.xlim(0,len(rollm))
plt.savefig(snakemake.output[0],bbox_inches='tight')
