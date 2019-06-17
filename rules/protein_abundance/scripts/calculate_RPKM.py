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
    sample2read_count[sample_id] = float(sample2read_count[sample_id])/1000000
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
        # divide by sequence length and number of million of reads
        counts = sample2sequence_accession2count[sample][sequence_accession]
        RPKM = counts/sample2read_count[sample]
        seq_length = record2aa_sequence_length[sequence_accession]
        norm_counts = counts/float(seq_length)
        #print("counts\t%s" % counts)
        #print("Seq length\t%s ==> %s" % (seq_length, norm_counts))
        reads_millions = sample2read_count[sample_id]
        #print("reads_millions\t%s" % reads_millions)
        #print(float(norm_counts), float(reads_millions))
        RPKM = float(norm_counts) / float(reads_millions)
        #print("RPKM\t%s" % RPKM)
        cursor.execute(sql_template, (sample, sequence_accession, n_hits, RPKM))
    conn.commit()

# index sequence accession
sql_index_1 = 'create index acc on sequence_counts (accession);'
cursor.execute(sql_index_1,)
conn.commit()
