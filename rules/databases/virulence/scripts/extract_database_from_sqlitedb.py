

import sqlite3
from Bio import SeqIO
conn = sqlite3.connect(snakemake.input[0])
cursor = conn.cursor()

record_dict = SeqIO.to_dict(SeqIO.parse(snakemake.input[1], "fasta"))

# get list of VFs from SwissProt and VFDB
sql = 'select distinct uniparc_accession from VF_table t1 inner join VF_databases t2 on t1.db_id=t2.db_id inner join uniparc2silix_90 t3 on t1.uniparc_id=t3.uniparc_id where db_name in ("VFDB", "SWISSPROT");'
uniparc_accession_list = [i[0] for i in cursor.execute(sql,).fetchall()]

records_subset = [record_dict[i] for i in uniparc_accession_list]

SeqIO.write(records_subset, snakemake.output[0], 'fasta')
