import pandas as pd
from Bio import SeqIO
fasta_db=snakemake.input[0]
handle=SeqIO.parse(fasta_db,'fasta')
desc=[]
gen_len=[]
ids=[]
for record in handle:
    desc.append(record.description.split(record.id+' ')[1].replace(' ','_').replace(',','-').replace('/','-').replace('[','_').replace(']','_').replace('{','_').replace('}','_'))
    gen_len.append(len(record.seq))
    acc=record.id.split('|')
    ids.append(acc[len(acc)-2].replace('.','_'))
genome_lengths_table=pd.DataFrame(gen_len,ids,columns=['genome_length'])
genome_lengths_table.to_csv(snakemake.output[0],sep='\t')
genome_names=pd.DataFrame(desc,ids,columns=['genome_names'])
genome_names.to_csv(snakemake.output[1],sep='\t')