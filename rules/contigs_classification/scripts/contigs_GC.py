
from Bio.SeqUtils import GC 
from Bio import SeqIO 

outfile = open(snakemake.output[0], 'w')
contig_records = SeqIO.parse(snakemake.input[0], "fasta")
for contig in contig_records:
    outfile.write("%s\t%s\n" % (contig.id, GC(contig.seq)))

