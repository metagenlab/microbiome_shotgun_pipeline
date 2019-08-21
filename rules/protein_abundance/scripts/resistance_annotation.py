
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
import pandas

conn = sqlite3.connect(snakemake.input["RPKM_database"])
cursor = conn.cursor()

protein_accession2aro_index = pandas.read_csv(snakemake.input["aro_index"], delimiter='\t').to_dict()["Protein Accession"]
aro_accession2description = pandas.read_csv(snakemake.input["aro"], delimiter='\t').to_dict()["Accession"]
protein_accession2aro_categories_index = pandas.read_csv(snakemake.input["aro_categories_index"], delimiter='\t').to_dict()["Protein Accession"]
RPKM_database = snakemake.input["RPKM_database"]

print(aro_accession2description)
