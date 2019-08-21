
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
import pandas

# ARO Accession
# Model ID
# Model Name
aro_index = pandas.read_csv(snakemake.input["aro_index"], delimiter='\t').set_index("Protein Accession")

# Name
# Description
aro = pandas.read_csv(snakemake.input["aro"], delimiter='\t').set_index("Accession")

# Resistance Mechanism
# AMR Gene Family
# Drug Class
aro_categories_index = pandas.read_csv(snakemake.input["aro_categories_index"], delimiter='\t').set_index("Protein Accession")

RPKM_database = snakemake.input["RPKM_database"]

conn = sqlite3.connect(snakemake.input["RPKM_database"])
cursor = conn.cursor()

# get uniparc accession list
sql = 'select distinct accession from sequence_counts'
cursor.execute(sql)
card_accession_list = [i[0] for i in cursor.fetchall()]

for accession in card_accession_list:
    aro_accession = aro_index.loc[accession, "ARO Accession"]
    aro_name = aro[aro_accession, "Name"]
    aro_description = aro[aro_accession, "Name"]
    resistance_mechanism = aro_categories_index[accession, "Resistance Mechanism"]
    AMR_family = aro_categories_index[accession, "AMR Gene Family"]
    drug_class = aro_categories_index[accession, "Drug Class"]
    print(accession)
    print(drug_class)
    print(resistance_mechanism)