import pandas

amrfinder_files = snakemake.input
sample_table = snakemake.params[0]

all_samples = pandas.read_csv(sample_table, sep="\t", index_col=0)

sample2group1 = []
sample2group2 = []

'''
0 Protein identifier                    NODE_10479_1	
1 Gene symbol	                        tet(M)	
2 Sequence name	                        tetracycline resistance ribosomal protection protein Tet(M)	
3 Scope	                                core
4 Element type	                        AMR	
5 Element subtype	                    AMR	
6 Class	                                TETRACYCLINE
7 Subclass	                            TETRACYCLINE
8 Method	                            PARTIALP
9 Target length	                        405
10 Reference sequence length            639	
11 % Coverage of reference sequence	    63.38
12 % Identity to reference sequence	    99.75
13 Alignment length	                    405
14 Accession of closest sequence	    WP_002414694.1
15 Name of closest sequence	            tetracycline resistance ribosomal protection protein Tet(M)
16 HMM id	                            NF012155.1
17 HMM description                      tetracycline resistance ribosomal protection protein Tet(M)
'''

header1 = ["tool",
           "sample", 
           "group_1",
           "group_2", 
           "contig_id",
           "contig_id_uniq",
           "ORF_id",
           "hit_accession", # data[14]
           "method", # data[8]
           "coverage_reference", # data[11]
           "identity", # data[12]
           "alignment_length", # data[13]
           "class", # data[6]
           "subclass", # data[7]
           "name", # data[1]
           "description", # data[2]
           "hmm_accession", # data[16]
           "hmm_description" # data[17]
           ]

header2 = ["aro", "drug_class"]
header3 = ["aro", "resistance_mechanism"]
header4 = ["aro", "gene_family"]

data = []

o_rgi_hits = open(snakemake.output[0], "w")
o_rgi_hits.write("\t".join(header1) + '\n')


for amrfinder_file in amrfinder_files:
    sample_name = amrfinder_file.split("/")[1]
    group_1 = all_samples.loc[sample_name, "group_1"]
    group_2 = all_samples.loc[sample_name, "group_2"]
    with open(amrfinder_file) as f:
        for n, row in enumerate(f):
            if n == 0:
                continue 
            else:
                data = row.rstrip().split("\t")
                # only keep AMR results
                if data[4] != "AMR":
                    continue
                if "plass" in rgi_file:
                    ## M03935:66:000000000-B656J:1:1101:14307:380 [Orf: 5677, 0, 135, -1, 1, 0] => renamed to ORF_X
                    tool = "plass"
                    contig_id = None
                    contig_id_uniq = None
                    orf_number = data[0]
                if "prodigal" in  rgi_file:
                    # ok 0 ORF_ID NODE_270_1 # 57 # 1982 # 1 # ID=270_1;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.398
                    tool = "prodigal"
                    contig_prefix, contig_number, orf_number = data[0].split(" ")[0].split("_")
                    contig_id = "%s_%s" % (contig_prefix, contig_number)
                    contig_id_uniq = "%s_%s" % (sample_name, contig_id)
                if "contig" in rgi_file:
                    # 
                hit_accession = data[14]
                method = data[8]
                coverage_reference = data[11]
                identity = data[12]
                alignment_length = data[13]
                amr_class = data[6]
                amr_subclass = data[7]
                name = data[1]
                description = data[2]
                hmm_accession = data[16]
                hmm_description = data[17]

                hit_accession = data[14]

                o_rgi_hits.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tool,
                                                                                                       sample_name,
                                                                                                       group_1,
                                                                                                       group_2,
                                                                                                       contig_id,
                                                                                                       contig_id_uniq,
                                                                                                       orf_number,
                                                                                                       hit_accession,
                                                                                                       method,
                                                                                                       coverage_reference,
                                                                                                       identity,
                                                                                                       alignment_length,
                                                                                                       amr_class,
                                                                                                       amr_subclass,
                                                                                                       name,
                                                                                                       description,
                                                                                                       hmm_accession,
                                                                                                       hmm_description
                                                                                                        ))


