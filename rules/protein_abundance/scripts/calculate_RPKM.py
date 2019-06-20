#!/usr/bin/env python

from Bio import SeqIO
import sqlite3

conn = sqlite3.connect(snakemake.output[0])
cursor = conn.cursor()

sample2read_count={}
# get total number of reads if R1 and R2
for R1_count in snakemake.input["R1_counts"]:
    with open(R1_count, 'r') as f:
        sample_id = R1_count.split("/")[1]
        sample2read_count[sample_id] = int(f.readline().rstrip().split("\t")[1])
for R2_count in snakemake.input["R2_counts"]:
    with open(R2_count, 'r') as f:
        sample_id = R2_count.split("/")[1]
        sample2read_count[sample_id] += int(f.readline().rstrip().split("\t")[1])
for sample_id in sample2read_count:
    print("reads", sample2read_count[sample_id])
    # get total number of reads
    sample2read_count[sample_id] = float(sample2read_count[sample_id])
print(sample2read_count)
# get sequence length for normalization by sequence length
records = SeqIO.parse(snakemake.input["reference_fasta"], 'fasta')
record2aa_sequence_length = {}
for record in records:
    record2aa_sequence_length[record.name] = len(record.seq)

# get number of match for each sequence
sample2sequence_accession2count = {}
for sample in snakemake.input["sample_list"]:
    sample_id = sample.split("/")[1]
    #print("sample_id", sample_id)
    sample2sequence_accession2count[sample_id] = {}
    with open(sample, 'r') as f:
        for line in f:
            data = line.rstrip().split("\t")
            if data[1] not in sample2sequence_accession2count[sample_id]:
                sample2sequence_accession2count[sample_id][data[1]] = 1
            else:
                sample2sequence_accession2count[sample_id][data[1]] += 1

sql = 'create table sequence_counts (sample varchar(200), accession varchar(200), n_hits INTEGER, RPKM INTEGER)'
cursor.execute(sql,)

# load data into database
sql_template = 'insert into sequence_counts values (?, ?, ?, ?)'
for sample in sample2sequence_accession2count:
    for sequence_accession in sample2sequence_accession2count[sample]:
        n_hits = sample2sequence_accession2count[sample][sequence_accession]

        # multiply by 1000000 to get reads per million
        # multiply by 1000 to calculate reads per kilobase
        nominat = sample2sequence_accession2count[sample][sequence_accession] * 1000000 * 1000

        # multiply protein length by 3 to get gene length
        seq_length_nucl = record2aa_sequence_length[sequence_accession]
        reads_millions = sample2read_count[sample_id]
        # multiply library size by gene length
        denominat = seq_length_nucl * reads_millions

        # calculate ratio
        RPKM = float(nominat)/float(denominat)

        cursor.execute(sql_template, (sample, sequence_accession, n_hits, RPKM))
    conn.commit()

# index sequence accession
sql_index_1 = 'create index acc on sequence_counts (accession);'
cursor.execute(sql_index_1,)
conn.commit()
