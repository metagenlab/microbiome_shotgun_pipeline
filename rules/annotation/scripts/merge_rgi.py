import pandas

rgi_files = snakemake.input["rgi_output"]
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
aro2accession = pandas.read_csv(card_annotation,
                                delimiter='\t',
                                names=["accession","accession_no_version", "ARO", "description"]).set_index("ARO").to_dict()["accession"]


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
           "contig_id_uniq"
           "ORF_id",
           "cut_Off",
           "pass_bitscore",
           "best_hit_bitscore",
           "best_hit_ARO",
           "best_identities",
           "ARO",
           "model_type",
           "percentage_length",
           "accession",
           "silix_98",
           "silix_98_decription",
           "silix_95",
           "silix_95_decription",
           "silix_90",
           "silix_90_decription",
           "silix_80",
           "silix_80_decription"]

header2 = ["aro", "drug_class"]
header3 = ["aro", "resistance_mechanism"]
header4 = ["aro", "gene_family"]

data = []

o_rgi_hits = open(snakemake.output[0], "w")
o_rgi_hits.write("\t".join(header1) + '\n')

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
                    ## M03935:66:000000000-B656J:1:1101:14307:380 [Orf: 5677, 0, 135, -1, 1, 0] => renamed ORF_X
                    tool = "plass"
                    contig_id = None
                    unique_id = None
                    orf_number = data[0] # .split("[Orf: ")[1].split(",")[0]
                if "prodigal" in  rgi_file:
                    # ok 0 ORF_ID NODE_270_1 # 57 # 1982 # 1 # ID=270_1;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.398
                    tool = "prodigal"
                    contig_prefix, contig_number, orf_number = data[0].split(" ")[0].split("_")
                    contig_id = "%s_%s" % (contig_prefix, contig_number)
                    contig_uniq = "%s_%s" % (sample_name, contig_id)
                if "contig" in rgi_file:
                    tool = "contigs"
                    if not "NODE_" in data[0]:
                        print("Skipping non protein homolog model data:", data[0])
                        continue
                    contig_prefix, contig_number, orf_number = data[0].split(" ")[0].split("_")
                    contig_id = "%s_%s" % (contig_prefix, contig_number)
                    contig_uniq = "%s_%s" % (sample_name, contig_id)
                                   
                cut_off = data[5]
                pass_bitscore = data[6]
                best_hit_bitscore = data[7]
                best_hit_ARO = data[8]
                best_identities = data[9]
                ARO = "ARO:%s" % data[10]
                model_type = data[11]
                percentage_length = data[20]
            
                if model_type != 'protein homolog model':
                    continue

                try:
                    accession = aro2accession[ARO]
                except KeyError:
                    print(data)
                    continue
                silix_98 = accession2silix_98[accession]
                silix_98_description = silix_98_2_annotation[silix_98]

                silix_95 = accession2silix_95[accession]
                silix_95_description = silix_95_2_annotation[silix_95]

                silix_90 = accession2silix_90[accession]
                silix_90_description = silix_90_2_annotation[silix_90]

                silix_80 = accession2silix_80[accession]
                silix_80_description = silix_80_2_annotation[silix_80]


                o_rgi_hits.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tool,
                                                                                              sample_name,
                                                                                               group_1,
                                                                                                group_2,
                                                                                                contig_id,
                                                                                                contig_uniq,
                                                                                                orf_number,
                                                                                                cut_off,
                                                                                                pass_bitscore,
                                                                                                best_hit_bitscore,
                                                                                                best_hit_ARO,
                                                                                                best_identities,
                                                                                                ARO,
                                                                                                model_type,
                                                                                                percentage_length,
                                                                                                accession,
                                                                                                silix_98,
                                                                                                silix_98_description,
                                                                                                silix_95,
                                                                                                silix_95_description,
                                                                                                silix_90,
                                                                                                silix_90_description,
                                                                                                silix_80,
                                                                                                silix_80_description
                                                                                           ))

