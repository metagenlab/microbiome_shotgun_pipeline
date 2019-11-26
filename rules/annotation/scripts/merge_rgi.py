import pandas

rgi_files = snakemake.input
sample_table = snakemake.params[0]

all_samples = pandas.read_csv(sample_table, sep="\t", index_col=0)

sample2group1 = []
sample2group2 = []

# ok 0 ORF_ID NODE_270_1 # 57 # 1982 # 1 # ID=270_1;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.398
# 1 Contig	
# 2 Start	
# 3 Stop	
# 4 Orientation	
# ok 5 Cut_Off	Strict
# ok 6 Pass_Bitscore 1200	
# ok 7 Best_Hit_Bitscore 1281.54	
# ok 8 Best_Hit_ARO tetQ	
# ok 9 Best_Identities	96.41
# ok 10 ARO 3000191	
# ok 11 Model_type	protein homolog model
# 12 SNPs_in_Best_Hit_ARO	macrolide antibiotic; fluoroquinolone antibiotic
# 13 Other_SNPs	
# ok 14 Drug Class	
# ok 15 Resistance Mechanism	antibiotic target alteration
# ok 16 AMR Gene Family	
# 17 Predicted_DNA	
# 18 Predicted_Protein	
# 19 CARD_Protein_Sequence	
# ok 20 Percentage Length of Reference Sequence	
# 21 ID	
# 22 Model_ID 1183	
# 23 Nudged	
# 24 Note


# single match only

# ok 0 ORF_ID NODE_270_1 # 57 # 1982 # 1 # ID=270_1;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.398

## PLASS: M03935:66:000000000-B656J:1:1101:14307:380 [Orf: 5677, 0, 135, -1, 1, 0]

# ok 5 Cut_Off	Strict
# ok 6 Pass_Bitscore 1200	
# ok 7 Best_Hit_Bitscore 1281.54	
# ok 8 Best_Hit_ARO tetQ	
# ok 9 Best_Identities	96.41
# ok 10 ARO 3000191	
# ok 11 Model_type	protein homolog model
# ok 20 Percentage Length of Reference Sequence	

## multiple match possibles ###

# ok 14 Drug Class	
# ok 15 Resistance Mechanism	antibiotic target alteration
# ok 16 AMR Gene Family	


header1 = ["tool",
           "sample", 
           "group_1",
           "group_2", 
           "contig_id",
           "ORF_id",
           "cut_Off",
           "pass_bitscore",
           "best_hit_bitscore",
           "best_hit_ARO",
           "best_identities",
           "ARO",
           "model_type",
           "percentage_length"]

header2 = ["aro", "drug_class"]
header3 = ["aro", "resistance_mechanism"]
header4 = ["aro", "gene_family"]

data = []

o_rgi_hits = open(snakemake.output[0], "w")
o_rgi_hits.write("\t".join(header1) + '\n')

o_aro2mechanism = open(snakemake.output[1], "w")
o_aro2mechanism.write("\t".join(header2) + '\n')

o_aro2gene_family = open(snakemake.output[2], "w")
o_aro2gene_family.write("\t".join(header3) + '\n')

o_aro2drug_class = open(snakemake.output[3], "w")
o_aro2drug_class.write("\t".join(header4) + '\n')

aro2mechanism = {}
aro2gene_family = {}
aro2drug_class = {}

for rgi_file in rgi_files:
    sample_name = rgi_file.split("/")[1]
    group_1 = all_samples.loc[sample_name, "group_1"]
    group_2 = all_samples.loc[sample_name, "group_2"]
    with open(rgi_file) as f:
        for n, row in enumerate(f):
            if n == 0:
                continue 
            else:
                
                data = row.rstrip().split("\t")
                if "plass" in rgi_file:
                    ## M03935:66:000000000-B656J:1:1101:14307:380 [Orf: 5677, 0, 135, -1, 1, 0]
                    tool = "plass"
                    contig_id = None
                    orf_number = data[0].split("[Orf: ")[1].split(",")[0]
                elif "prodigal" in  rgi_file:
                    # ok 0 ORF_ID NODE_270_1 # 57 # 1982 # 1 # ID=270_1;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.398
                    tool = "prodigal"
                    contig_prefix, contig_number, orf_number = data[0].split(" ")[0].split("_")
                    contig_id = "%s_%s" % (contig_prefix, contig_number)

                cut_off = data[5]
                pass_bitscore = data[6]
                best_hit_bitscore = data[7]
                best_hit_ARO = data[8]
                best_identities = data[9]
                ARO = "ARO:%s" % data[10]
                model_type = data[11]
                percentage_length = data[20]

                o_rgi_hits.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tool,
                                                                                           sample_name,
                                                                                           group_1,
                                                                                           group_2,
                                                                                           contig_id,
                                                                                           orf_number,
                                                                                           cut_off,
                                                                                           pass_bitscore,
                                                                                           best_hit_bitscore,
                                                                                           best_hit_ARO,
                                                                                           best_identities,
                                                                                           ARO,
                                                                                           model_type,
                                                                                           percentage_length
                                                                                           ))

                if ARO not in aro2mechanism:
                    drug_class_list = data[14].split(";")
                    mechanism_list = data[15].split(";")
                    gene_family_list = data[16].split(";")

                    aro2mechanism[ARO] = mechanism_list
                    aro2gene_family[ARO] = gene_family_list
                    aro2drug_class[ARO] = drug_class_list
for aro in aro2mechanism:
    for mechanism in aro2mechanism[aro]:
        o_aro2mechanism.write("%s\t%s\n" % (aro, mechanism))

    for family in aro2gene_family[aro]:
        o_aro2gene_family.write("%s\t%s\n" % (aro, family))   

    for drug_class in aro2drug_class[aro]:
        o_aro2drug_class.write("%s\t%s\n" % (aro, drug_class))