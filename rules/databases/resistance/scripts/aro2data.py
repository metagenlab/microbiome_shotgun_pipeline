from Bio import SeqIO 

o = open(snakemake.output[0], "w") # accession2AMR_family.tab
p = open(snakemake.output[1], "w") # accession2drug_class.tab
q = open(snakemake.output[2], "w") # accession2resistance_mechanism.tab

aro_index = snakemake.input[0]
aro_categories_index = snakemake.input[1]
CARD_homolog_model_fasta = snakemake.input[2]


def get_accession2aro(aro_index):

    protein_accession2ARO = {}
    with open(aro_index, 'r') as f:
        for n, row in enumerate(f):
            if n == 0:
                continue
            else:
                data = row.rstrip().split("\t")
                protein_accession2ARO[data[6]] = data[0]
    return protein_accession2ARO


def get_aro2category_data(aro_categories_index, protein_accession2ARO):
    aro2category_data = {}
    with open(aro_categories_index, 'r') as f:
        for	n, row in enumerate(f):
            if n ==	0:
                continue
            else:
                data = row.rstrip().split("\t")
                protein_accession = data[0]
                if protein_accession == 'N/A':
                    print("Skipping mutation conferring resistance:", data[2])
                    continue
                ARO = protein_accession2ARO[protein_accession]
                aro2category_data[ARO] = {}

                AMR_family_list = data[2].split(";")
                drug_class_list = data[3].split(";")
                resistance_mechanism_list = data[4].split(";")

                aro2category_data[ARO]["AMR_family_list"] = AMR_family_list
                aro2category_data[ARO]["drug_class_list"] = drug_class_list
                aro2category_data[ARO]["resistance_mechanism_list"] = resistance_mechanism_list
    return aro2category_data


accession2aro = get_accession2aro(aro_index)
aro2category_data = get_aro2category_data(aro_categories_index, accession2aro)

records = SeqIO.parse(CARD_homolog_model_fasta, 'fasta')

for record in records:
    protein_accession = record.id
    protein_accession_no_version = protein_accession.split(".")[0]
    ARO = accession2aro[protein_accession]
    
    AMR_family_list = aro2category_data[ARO]["AMR_family_list"]
    drug_class_list = aro2category_data[ARO]["drug_class_list"]
    resistance_mechanism_list = aro2category_data[ARO]["resistance_mechanism_list"]

    for family in AMR_family_list:
        o.write("%s\t%s\t%s\t%s\n" % (protein_accession, protein_accession_no_version, ARO, family))
    for drug_class in drug_class_list:
        p.write("%s\t%s\t%s\t%s\n" % (protein_accession, protein_accession_no_version, ARO, drug_class))
    for mechanism in resistance_mechanism_list:
        q.write("%s\t%s\t%s\t%s\n" % (protein_accession, protein_accession_no_version, ARO, mechanism))