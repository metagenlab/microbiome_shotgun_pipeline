
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
import pandas

# output

o = open(snakemake.output[0], 'w')

RPKM_table = snakemake.input["RPKM_table"]

card_annotation = snakemake.input["card_annotation"]
silix_98 = snakemake.input["silix_98"]
silix_95 = snakemake.input["silix_95"]
silix_90 = snakemake.input["silix_90"]
silix_80 = snakemake.input["silix_80"]
silix_98_annotation = snakemake.input["silix_98_annotation"]
silix_95_annotation = snakemake.input["silix_95_annotation"]
silix_90_annotation = snakemake.input["silix_90_annotation"]
silix_80_annotation = snakemake.input["silix_80_annotation"]

accession2silix_98 = pandas.read_csv(silix_98,
                                     delimiter='\t',
                                     names=["silix_name","accession"]).set_index("accession").to_dict()["silix_name"]

accession2silix_95 = pandas.read_csv(silix_95,
                                     delimiter='\t',
                                     names=["silix_name","accession"]).set_index("accession").to_dict()["silix_name"]

accession2silix_90 = pandas.read_csv(silix_90,
                                     delimiter='\t',
                                     names=["silix_name","accession"]).set_index("accession").to_dict()["silix_name"]

accession2silix_80 = pandas.read_csv(silix_80,
                                     delimiter='\t',
                                     names=["silix_name","accession"]).set_index("accession").to_dict()["silix_name"]

silix_98_2_annotation = pandas.read_csv(silix_98_annotation,
                                     delimiter='\t',
                                     names=["silix_name","description"]).set_index("silix_name").to_dict()["description"]

silix_95_2_annotation = pandas.read_csv(silix_95_annotation,
                                     delimiter='\t',
                                     names=["silix_name","description"]).set_index("silix_name").to_dict()["description"]

silix_90_2_annotation = pandas.read_csv(silix_90_annotation,
                                     delimiter='\t',
                                     names=["silix_name","description"]).set_index("silix_name").to_dict()["description"]

silix_80_2_annotation = pandas.read_csv(silix_80_annotation,
                                     delimiter='\t',
                                     names=["silix_name","description"]).set_index("silix_name").to_dict()["description"]

# ACT97415.1	ACT97415	ARO:3002999	CblA beta-lactamase
accession2aro = pandas.read_csv(card_annotation,
                                delimiter='\t',
                                names=["accession","accession_no_version", "ARO", "description"]).set_index("accession").to_dict()["ARO"]


with open(RPKM_table, 'r') as f:
    for n, row in enumerate(f):
        data = row.split("\t")
        if n == 0:
            # header 
            # sample sequence_accession n_hits RPKM group_1 group_2
            o.write(row.rstrip() + '\tARO\tsilix_98\tsilix_98_description\tsilix_95\tsilix_95_description\tsilix_90\tsilix_90_description\tsilix_80\tsilix_80_description\n')
        else:
            accession = data[1]
            aro = accession2aro[accession]

            silix_98 = accession2silix_98[accession]
            silix_98_description = silix_98_2_annotation[silix_98]

            silix_95 = accession2silix_95[accession]
            silix_95_description = silix_95_2_annotation[silix_95]

            silix_90 = accession2silix_90[accession]
            silix_90_description = silix_90_2_annotation[silix_90]

            silix_80 = accession2silix_80[accession]
            silix_80_description = silix_80_2_annotation[silix_80]
            o.write(row.rstrip()+ f'\t{aro}\t{silix_98}\t{silix_98_description}\t{silix_95}\t{silix_95_description}\t{silix_90}\t{silix_90_description}\t{silix_80}\t{silix_80_description}\n')

