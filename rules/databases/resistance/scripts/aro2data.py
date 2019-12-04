o = open(snakemake.output[0], "w") # accession2AMR_family.tab
p = open(snakemake.output[1], "w") # accession2drug_class.tab
q = open(snakemake.output[2], "w") # accession2resistance_mechanism.tab

aro_index = snakemake.input[0]
aro_categories_index = snakemake.input[1]

protein_accession2ARO = {}
with open(aro_index, 'r') as f:
    for n, row in enumerate(f):
        if n == 0:
            continue
        else:
            data = row.rstrip().split("\t")
            protein_accession2ARO[data[6]] = data[0]

with open(aro_categories_index, 'r') as f:
    for	n, row in enumerate(f):
        if n ==	0:
            continue
        else:
            data = row.rstrip().split("\t")
            protein_accession = data[0]
            protein_accession_no_version = data[0].split(".")[0]
            ARO = protein_accession2ARO[protein_accession]
            AMR_family_list = data[2].split(";")
            for family in AMR_family_list:
                o.write("%s\t%s\t%s\t%s\n" % (protein_accession, protein_accession_no_version, ARO, family))
            drug_class_list = data[3].split(";")
            for drug_class in drug_class_list:
                p.write("%s\t%s\t%s\t%s\n" % (protein_accession, protein_accession_no_version, ARO, drug_class))
            resistance_mechanism_list = data[4].split(";")
            for mechanism in resistance_mechanism_list:
                q.write("%s\t%s\t%s\t%s\n" % (protein_accession, protein_accession_no_version, ARO, mechanism))