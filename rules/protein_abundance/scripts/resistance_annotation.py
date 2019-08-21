
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
import pandas

conn = sqlite3.connect(snakemake.input["RPKM_database"])
cursor = conn.cursor()

aro_index = pandas.read_csv(snakemake.input["aro_index"], delimiter='\t').set_index("Protein Accession")
aro = pandas.read_csv(snakemake.input["aro"], delimiter='\t').set_index("Accession")
aro_categories_index = pandas.read_csv(snakemake.input["aro_categories_index"], delimiter='\t').set_index("Protein Accession")
RPKM_database = snakemake.input["RPKM_database"]

print(aro)
