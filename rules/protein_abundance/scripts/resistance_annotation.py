
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
import pandas

# output

o = open(snakemake.output[0], 'w')

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

sql1 = 'create table accession2aro (protein_accession varchar(200), aro_accession varchar(200), aro_name varchar(200), aro_description TEXT, resistance_mechanism TEXT, AMR_family TEXT)'
cursor.execute(sql1)
sql2 = 'create table aro_accession2drug_class (aro_accession varchar(200), drug_class TEXT)'
cursor.execute(sql2)

aro_accession2drug_class_list = {}

sql_template = 'insert into accession2aro values (?, ?, ?, ?, ?, ?)'
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
    except:
        # some missing sequences in aro_categories_index
        # use another seq to retrieve aro annotation
        alternative_sequences = aro_index.loc[(aro_index["ARO Accession"] == aro_accession) & (aro_index.index != accession)].index.tolist()
        resistance_mechanism = aro_categories_index.loc[alternative_sequences[0], "Resistance Mechanism"]
        AMR_family = aro_categories_index.loc[alternative_sequences[0], "AMR Gene Family"]
        drug_class = aro_categories_index.loc[alternative_sequences[0], "Drug Class"]
    if aro_accession not in aro_accession2drug_class_list:
        aro_accession2drug_class_list[aro_accession] = drug_class.split(";")
    
    cursor.execute(sql_template, (accession,
                                  aro_accession,
                                  aro_name,
                                  aro_description,
                                  resistance_mechanism,
                                  AMR_family))

sql_template_2 = 'insert into aro_accession2drug_class values (?, ?)'

for aro in aro_accession2drug_class_list:
    for drug_class in aro_accession2drug_class_list[aro]:
        cursor.execute(sql_template_2, (aro, drug_class))

conn.commit()   
o.write("ok")
o.close()