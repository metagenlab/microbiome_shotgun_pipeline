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

sql = 'create table uniparc_accession2annotation (accession varchar(200), silix_accession varchar(200), annotation TEXT)'

sql_taxonomy = 'insert into uniparc_accession2genus values (?, ?)'

cursor.execute(sql,)

# load data into database
sql_template = 'insert into uniparc_accession2annotation values (?, ?, ?)'
acc2genus_ok = []

# get uniparc accession list
sql = 'select distinct accession from sequence_counts'
cursor.execute(sql)
uniparc_accession_list = [i[0] for i in cursor.fetchall()]

for uniparc_accession in uniparc_accession_list:

    silix_acc = uniparc_accession2silix_acc_and_consensus_annotation[uniparc_accession][0]
    description = uniparc_accession2silix_acc_and_consensus_annotation[uniparc_accession][1]
    
    cursor.execute(sql_template, (uniparc_accession,
                                    silix_acc,
                                    description))

    if uniparc_accession not in acc2genus_ok:
        genus_list = uniparc_accession2genus_list[uniparc_accession]
        for genus in genus_list:
            cursor.execute(sql_taxonomy, [uniparc_accession, genus])
        acc2genus_ok.append(uniparc_accession)
conn.commit()

# index sequence accession
sql_index_1 = 'create index acc3 on uniparc_accession2annotation (accession);'
cursor.execute(sql_index_1,)

sql_index_2 = 'create index acc2 on uniparc_accession2genus (accession);'
cursor.execute(sql_index_2,)

conn.commit()
