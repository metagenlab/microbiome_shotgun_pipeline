#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
import pandas
import re

conn = sqlite3.connect(snakemake.output[0])
cursor = conn.cursor()

sample_table = snakemake.params[0]

all_samples = pandas.read_csv(sample_table, sep="\t", index_col=0)

sample2read_count={}
# get total number of reads if R1 and R2
for flash_log in snakemake.input["flash_log"]:
    sample_id = flash_log.split("/")[1]
    with open(flash_log, 'r') as f:
        for line in f:
            if "Total pairs" in line:
                total_reads = re.findall('\d+', line)[0]
                # [FLASH]     Total pairs:      4704430
                sample2read_count[sample_id] = float(total_reads)
# for R2_count in snakemake.input["R2_counts"]:
#   with open(R2_count, 'r') as f:
#        sample_id = R2_count.split("/")[1]
#        sample2read_count[sample_id] += int(f.readline().rstrip().split("\t")[1])

print(sample2read_count)
# get sequence length for normalization by sequence length
records = SeqIO.parse(snakemake.input["reference_fasta"], 'fasta')
record2aa_sequence_length = {}
for record in records:
    if '|' in record.id:
        acc = record.id.split("|")[1]
    else:
        acc = record.id
    record2aa_sequence_length[acc] = len(record.seq)
print(record2aa_sequence_length)
# get number of match for each sequence
sample2sequence_accession2count = {}
for sample in snakemake.input["sample_list"]:
    sample_id = sample.split("/")[1]
    #print("sample_id", sample_id)
    sample2sequence_accession2count[sample_id] = {}
    with open(sample, 'r') as f:
        for line in f:
            data = line.rstrip().split("\t")
            # gb|AAG07064.1|ARO:3003693|mexK
            # AAG07064.1
            if '|' in data[1]:
                hit_accession = data[1].split("|")[1]
            else:
                hit_accession = data[1]
            if hit_accession not in sample2sequence_accession2count[sample_id]:
                sample2sequence_accession2count[sample_id][hit_accession] = 1
            else:
                sample2sequence_accession2count[sample_id][hit_accession] += 1

sql = 'create table sequence_counts (sample varchar(300), accession varchar(300), n_hits INTEGER, RPKM FLOAT, group_1 varchar(300), group_2 varchar (300))'
cursor.execute(sql)

# load data into database
sql_template = 'insert into sequence_counts values (?, ?, ?, ?, ?, ?)'
acc2genus_ok = []
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
        # print("sample", sample)
        group_1 = all_samples.loc[sample, "group_1"]
        group_2 = all_samples.loc[sample, "group_2"]
        # print("RPKM", RPKM)
        # print("groups", group_1, group_2)
        cursor.execute(sql_template, (sample,
                                      sequence_accession,
                                      n_hits,
                                      RPKM,
                                      group_1,
                                      group_2))
    conn.commit()

# index sequence accession
sql_index_1 = 'create index acc on sequence_counts (accession);'
cursor.execute(sql_index_1,)

conn.commit()
