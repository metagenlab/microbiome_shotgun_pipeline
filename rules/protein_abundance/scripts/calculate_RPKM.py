#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
import pandas

conn = sqlite3.connect(snakemake.output[0])
cursor = conn.cursor()

sql = 'attach "%s" as virulencedb' % snakemake.input["virulence_database"]
cursor.execute(sql)

sql = 'create table uniparc_accession2genus (accession varchar(200), genus varchar(400))'
cursor.execute(sql)

sample_table = snakemake.params[0]

uniparc_accession2silix_acc_and_consensus_annotation = {}
sql = 'select uniparc_accession,silix_name,description from uniparc2silix_90 t1 inner join uniparc_consensus_annotation t2 on t1.uniparc_id=t2.uniparc_id;'
for row in cursor.execute(sql).fetchall():
    row = list(row)
    uniparc_accession2silix_acc_and_consensus_annotation[row[0]] = row[1:]

uniparc_accession2genus_list = {}
sql = 'select distinct uniparc_accession,genus_name from uniparc2silix_90 t1 inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id;'
for row in cursor.execute(sql).fetchall():
    row = list(row)
    if row[0] not in uniparc_accession2genus_list:
        uniparc_accession2genus_list[row[0]] = [row[1]]
    else:
        uniparc_accession2genus_list[row[0]].append(row[1])


all_samples = pandas.read_csv(sample_table, sep="\t", index_col=0)

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
            if data[1] not in sample2sequence_accession2count[sample_id]:
                sample2sequence_accession2count[sample_id][data[1]] = 1
            else:
                sample2sequence_accession2count[sample_id][data[1]] += 1

sql = 'create table sequence_counts (sample varchar(200), accession varchar(200), n_hits INTEGER, RPKM INTEGER, group_1 varchar(200), group_2 varchar (200), silix_accession varchar(200), annotation TEXT)'

sql_taxonomy = 'insert into uniparc_accession2genus values (?, ?)'
cursor.execute(sql,)

# load data into database
sql_template = 'insert into sequence_counts values (?, ?, ?, ?, ?, ?, ?, ?)'
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

        silix_acc = uniparc_accession2silix_acc_and_consensus_annotation[sequence_accession][0]
        description = uniparc_accession2silix_acc_and_consensus_annotation[sequence_accession][1]
        cursor.execute(sql_template, (sample,
                                      sequence_accession,
                                      n_hits,
                                      RPKM,
                                      all_samples.loc[sample, "group_1"],
                                      all_samples.loc[sample, "group_2"],
                                      silix_acc,
                                      description))
        if sequence_accession not in acc2genus_ok:
            genus_list = uniparc_accession2genus_list[sequence_accession]
            for genus in genus_list:
                cursor.execute(sql_taxonomy, [sequence_accession, genus])
            acc2genus_ok.append(sequence_accession)
    conn.commit()

# index sequence accession
sql_index_1 = 'create index acc on sequence_counts (accession);'
cursor.execute(sql_index_1,)

sql_index_2 = 'create index acc2 on uniparc_accession2genus (accession);'
cursor.execute(sql_index_2,)

conn.commit()
