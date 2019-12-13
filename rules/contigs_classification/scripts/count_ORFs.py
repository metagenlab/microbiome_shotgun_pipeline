
from Bio import SeqIO

prodigal_faa = snakemake.input[0]

o = open(snakemake.output[0], "w")

contig2count = {}

sample_name = prodigal_faa.split("/")[1]

for record in SeqIO.parse(prodigal_faa, "fasta"):
    contig_prefix, contig_number, orf_number = record.id.split(" ")[0].split("_")
    contig_id = f"{contig_prefix}_{contig_number}"
    if contig_id not in contig2count:
        contig2count[contig_id] = 1
    else:
        contig2count[contig_id] += 1

for contig in contig2count:
    o.write(f"%s\t%s\n" % (contig,
                           contig2count[contig]))