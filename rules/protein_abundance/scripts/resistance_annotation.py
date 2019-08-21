
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

print(aro)
problems = 0
for accession in card_accession_list:
    aro_accession = aro_index.loc[accession, "ARO Accession"]
    aro_name = aro.loc[aro_accession, "Name"]
    aro_description = aro.loc[aro_accession, "Description"]
    try:
        resistance_mechanism = aro_categories_index.loc[accession, "Resistance Mechanism"]
        AMR_family = aro_categories_index.loc[accession, "AMR Gene Family"]
        drug_class = aro_categories_index.loc[accession, "Drug Class"]
        #print(accession)
        # macrolide antibiotic;fluoroquinolone antibiotic;aminoglycoside antibiotic;lincosamide antibiotic;carbapenem;fosfomycin;cephalosporin;glycylcycline;bicyclomycin;penam;nucleoside antibiotic;tetracycline antibiotic;peptide antibiotic;acridine dye;oxazolidinone antibiotic;rifamycin antibiotic;diaminopyrimidine antibiotic;phenicol antibiotic;isoniazid;penem;benzalkonium chloride;rhodamine;antibacterial free fatty acids;nitroimidazole antibiotic
        #print(drug_class)
        #print(resistance_mechanism)
    exept: 
        print("problem with", aro_accession)
        problems+=1
print(problems)
