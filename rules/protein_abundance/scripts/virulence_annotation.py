#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
import pandas

conn = sqlite3.connect(snakemake.input["RPKM_database"])
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

sql = 'create table uniparc_accession2annotation (accession varchar(200), silix_accession varchar(200), annotation TEXT)'

sql_taxonomy = 'insert into uniparc_accession2genus values (?, ?)'

cursor.execute(sql,)

# load data into database
sql_template = 'insert into uniparc_accession2annotation values (?, ?, ?)'
acc2genus_ok = []
for sample in sample2sequence_accession2count:

        silix_acc = uniparc_accession2silix_acc_and_consensus_annotation[sequence_accession][0]
        description = uniparc_accession2silix_acc_and_consensus_annotation[sequence_accession][1]
        
        cursor.execute(sql_template, (sequence_accession,
                                      silix_acc,
                                      description))

        if sequence_accession not in acc2genus_ok:
            genus_list = uniparc_accession2genus_list[sequence_accession]
            for genus in genus_list:
                cursor.execute(sql_taxonomy, [sequence_accession, genus])
            acc2genus_ok.append(sequence_accession)
    conn.commit()

# index sequence accession
sql_index_1 = 'create index acc3 on uniparc_accession2annotation (accession);'
cursor.execute(sql_index_1,)

sql_index_2 = 'create index acc2 on uniparc_accession2genus (accession);'
cursor.execute(sql_index_2,)

conn.commit()
