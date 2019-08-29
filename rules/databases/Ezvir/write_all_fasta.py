db_path_ezvir=config["ezvir_database"]
from Bio import SeqIO
import os
handle=SeqIO.parse(snakemake.input[0],'fasta')
os.mkdir(os.path.join(db_path_ezvir,"index"))
for record in handle:
    acc=record.id.split('|')
    ids=acc[len(acc)-2].replace('.','_')
    genome_path=os.path.join(db_path_ezvir,"index",ids+".fa")
    SeqIO.write(record,genome_path,'fasta')