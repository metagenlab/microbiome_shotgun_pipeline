

aro2family_file = snakemake.input[0]
silix_file = snakemake.input[1]


def parse_silix(silix_node_file):
    silix_name2members = {}

    with open(silix_node_file, "r") as f:
        for row in f:
            data = row.rstrip().split("\t")
            silix_group = data[0]
            accession = data[1]
            if silix_group not in silix_name2members:
                silix_name2members[silix_group] = [accession]
            else:
                silix_name2members[silix_group].append(accession)

    return silix_name2members

def parse_silix_name2members(aro2family_file):
    prot_accession2family = {}

    with open(aro2family_file, "r") as f:
        for row in f:
            # AAN84550.1	AAN84550	ARO:3002276	VIM beta-lactamase
            data = row.rstrip().split("\t")
            protein_acc = data[0]
            protein_acc_no_version = data[1]
            aro = data[2]
            family = data[3]
            prot_accession2family[protein_acc] = family

    return prot_accession2family


def get_consensus_annotation(silix_name2members, prot_accession2family):
    silix_name2consensus_family = {}

    for silix_group in silix_name2members:
        acc_list = silix_name2members[silix_group]
        fam_list = list(set([prot_accession2family[i] for i in acc_list]))
        silix_name2consensus_family[silix_group] = ';'.join(fam_list)

    return silix_name2consensus_family


silix_name2members = parse_silix(silix_file)
prot_accession2family = parse_silix_name2members(aro2family_file)

silix_name2consensus_family = get_consensus_annotation(silix_name2members, prot_accession2family)

with open(snakemake.output[0], 'w') as f:
    for silix_group in silix_name2consensus_family:
        f.write("%s\t%s\n" % (silix_group, silix_name2consensus_family[silix_group]))