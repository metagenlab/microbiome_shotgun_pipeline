#!/usr/bin/env python

from Bio import SearchIO
import re 
import sqlite3
import pandas

rgi = snakemake.input["rgi"]
core_genes = snakemake.input["core_genes"]
plsdb = snakemake.input["plsdb"]
COG_mobilome = snakemake.input["COG_mobilome"]
deepvirfinder = snakemake.input["deepvirfinder"]
depth = snakemake.input["depth"]
GC = snakemake.input["GC"]
n_ORFs = snakemake.input["n_ORFs"]

plsdb_info_table = snakemake.params[0]
all_samples = pandas.read_csv(snakemake.params[1], sep="\t", index_col=0)

conn = sqlite3.connect(snakemake.output[0])
cursor = conn.cursor()

class Contig:

    def __init__(self,
                 rgi, 
                 core_genes, 
                 plsdb, 
                 COG_mobilome, 
                 deepvirfinder, 
                 depth, 
                 GC,
                 n_ORFs,
                 plsdb_info_table):

        self.rgi = rgi # ok 
        self.core_genes = core_genes # ok
        self.plsdb = plsdb
        self.COG_mobilome = COG_mobilome # ok
        self.deepvirfinder = deepvirfinder # ok
        self.depth = depth # ok
        self.GC = GC # ok
        self.n_ORFs = n_ORFs
        self.plsdb_info_table = plsdb_info_table
        self.plsdb_entry2length = {} 
        self.sample2plsdb = {} 
       # parse plasmid table to get plasmid length 
        self.parse_plsdb_info_table()

        self.sample2contig = {}
        self.sample2core_genes = {} 
        self.sample2COG_mobilome = {} 
        self.sample2rgi = {} 
        self.aro2drug_class_list = {}
        self.aro2resistance_mechanism_list = {}

    def parse_depth(self):
        # contigName	contigLen	totalAvgDepth	Amsterdam_RAW.bam	Amsterdam_RAW.bam-var
        # samples/Amsterdam_RAW/bwa/Amsterdam_RAW_assembled/Amsterdam_RAW_metabat2_depth.txt
        for depth_file in self.depth:
            sample = depth_file.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
            with open(depth_file) as f:
                for n, row in enumerate(f):
                    if n == 0:
                        continue
                    data = row.rstrip().split("\t")
                    contig_name = data[0]
                    contig_length = data[1]
                    average_depth = data[3]
                    var_depth = data[4]
                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    self.sample2contig[sample][contig_name]["contig_length"] = contig_length
                    self.sample2contig[sample][contig_name]["average_depth"] = average_depth
                    self.sample2contig[sample][contig_name]["var_depth"] = var_depth

    def parse_GC(self,):
        for gc_file in self.GC:
            sample = gc_file.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
            with open(gc_file) as f:
                for n, row in enumerate(f):
                    data = row.rstrip().split("\t")
                    contig_name = data[0]
                    gc = data[1] 
                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    self.sample2contig[sample][contig_name]["GC"] = gc


    def parse_n_ORFs(self,):
        for ORF_file in self.n_ORFs:
            sample = ORF_file.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
            with open(ORF_file) as f:
                for n, row in enumerate(f):
                    data = row.rstrip().split("\t")
                    contig_name = data[0]
                    n_ORFs = data[1] 
                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    self.sample2contig[sample][contig_name]["n_ORFs"] = n_ORFs


    def parse_deepvirfinder(self,):
        # name	len	score	pvalue
        for deepvirfinder_file in self.deepvirfinder:
            sample = deepvirfinder_file.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
            with open(deepvirfinder_file) as f:
                for n, row in enumerate(f):
                    if n == 0:
                        continue
                    data = row.rstrip().split("\t")
                    contig_name = data[0]
                    deepvirfinder_score = data[2]
                    deepvirfinder_pval = data[3]
                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    self.sample2contig[sample][contig_name]["deepvirfinder"] = {"score" : deepvirfinder_score, "pvalue" : deepvirfinder_pval}
    

    def parse_cog_mobilome(self,):
        '''

        0    qseqid     NODE_4_14	
        1    sseqid     CDD:225484	
        2    pident     23.121	
        3    length     173	
        4    mismatch   118	
        5    gapopen    5	
        6    qstart     59	
        7    qend       217	
        8    sstart     42	
        9    send       213	
        10   evalue     1.55e-23	
        11   bitscore   86.7
        '''

        contig_regex = re.compile("(NODE_[0-9]+).*")

        for COG_mobilome in self.COG_mobilome:
            sample = COG_mobilome.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
            if sample not in self.sample2COG_mobilome:
                self.sample2COG_mobilome[sample] = {}  
            with open(COG_mobilome) as f:
                for n, row in enumerate(f):
                    data = row.rstrip().split("\t")
                    ORF_id = data[0]
                    COG_cdd_id = data[1]
                    pident = data[2]
                    length = data[3]
                    mismatch = data[4]
                    gapopen = data[5]
                    qstart = data[6]
                    qend = data[7]
                    sstart = data[8]
                    send = data[9]
                    evalue = data[10]
                    bitscore = data[11]

                    contig_name = contig_regex.search(ORF_id).group(1)

                    # only keep first hit, else skip
                    if ORF_id not in self.sample2COG_mobilome:
                        self.sample2COG_mobilome[sample][ORF_id] = data 
                    else:
                        continue
                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    if "COG_mobilome" not in self.sample2contig[sample][contig_name]:
                        self.sample2contig[sample][contig_name]["COG_mobilome"] = 1
                    else:
                        self.sample2contig[sample][contig_name]["COG_mobilome"] += 1
 

    def parse_core_genes(self,):
        contig_regex = re.compile("(NODE_[0-9]+)_.*")

        for core_genes in self.core_genes:
            sample = core_genes.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
            if sample not in self.sample2core_genes:
                self.sample2core_genes[sample] = []  
            records = SearchIO.parse(core_genes, 'hmmer3-tab')
            for record in records:
                for hit in record.hits:
                    ORF_id = hit.id
                    HMM_id = record.id
                    evalue = hit.evalue
                    contig_name = contig_regex.search(ORF_id).group(1)

                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    # count number of core genes for each contig 
                    if "core_genes" not in self.sample2contig[sample][contig_name]:
                        self.sample2contig[sample][contig_name]["core_genes"] = 1
                    else:
                        self.sample2contig[sample][contig_name]["core_genes"] += 1
                    
                    # store detailed results 
                    self.sample2core_genes[sample].append([contig_name, ORF_id, HMM_id, evalue]) 

    def parse_rgi(self,):

        '''
        * 0	ORF_ID	NODE_44_2 # 428 # 2401 # -1 # ID=44_2;partial=00;start_type=GTG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;gc_cont=0.397
        1	Contig	
        2	Start	
        3	Stop	
        4	Orientation	
        *5	Cut_Off	Strict
        *6	Pass_Bitscore	1200
        *7	Best_Hit_Bitscore	1314.29
        *8	Best_Hit_ARO	tetQ
        *9	Best_Identities	96.35
        *10	ARO	3000191
        *11	Model_type	protein homolog model
        12	SNPs_in_Best_Hit_ARO	n/a
        13	Other_SNPs	n/a
        *14	Drug Class	macrolide antibiotic; lincosamide antibiotic; streptogramin antibiotic
        *15	Resistance Mechanism	antibiotic target alteration; antibiotic target replacement
        *16	AMR Gene Family	tetracycline-resistant ribosomal protection protein
        17	Predicted_DNA	
        18	Predicted_Protein	-
        19	CARD_Protein_Sequence	-
        *20	Percentage Length of Reference Sequence	100
        21	ID	gnl|BL_ORD_ID|1116|hsp_num:0
        22	Model_ID	1183
        23	Nudged	-
        24	Note	-
         '''

        contig_regex = re.compile("(NODE_[0-9]+)_.*")

        for rgi in self.rgi:
            sample = rgi.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
            if sample not in self.sample2rgi:
                self.sample2rgi[sample] = [] 
            with open(rgi) as f:
                for n, row in enumerate(f):
                    if n == 0:
                        continue
                    data = row.rstrip().split("\t")
                    # print(data)
                    ORF_id = data[0].split("#")[0] 
                    contig_name = contig_regex.search(ORF_id).group(1)
                    cut_off = data[5] 
                    Pass_Bitscore =  data[6]
                    Best_Hit_Bitscore = data[7]
                    Best_Hit_ARO =  data[8]
                    Best_Identities =  data[9]
                    ARO = data[10]
                    Model_type = data[11]
                    drug_class_list = data[14].split(';') 
                    resistance_mechanism_list = data[15].split(';') 
                    AMR_gene_family =  data[16]
                    percentage_length_of_reference_sequence =  data[20]

                    if ARO not in  self.aro2drug_class_list:
                         self.aro2drug_class_list[ARO] =  drug_class_list
                    if ARO not in  self.aro2resistance_mechanism_list:
                         self.aro2resistance_mechanism_list[ARO] =  resistance_mechanism_list

                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    
                   # only consider homolog model 
                    if Model_type == 'protein homolog model':
                        if "rgi" not in self.sample2contig[sample][contig_name]:
                            self.sample2contig[sample][contig_name]["rgi"] = 1
                        else:
                            self.sample2contig[sample][contig_name]["rgi"] += 1
                    
                        self.sample2rgi[sample].append([contig_name,
                                                        ORF_id,
                                                        cut_off,
                                                        Pass_Bitscore,
                                                        Best_Hit_Bitscore,
                                                        Best_Hit_ARO,
                                                        Best_Identities,
                                                        ARO,
                                                        Model_type,
                                                        AMR_gene_family,
                                                        percentage_length_of_reference_sequence])



    def parse_plsdb_info_table(self,):
        '''
        0	UID_NUCCORE
        1	ACC_NUCCORE
        2	Description_NUCCORE
        3	CreateDate_NUCCORE
        4	Topology_NUCCORE
        5	Completeness_NUCCORE
        6	TaxonID_NUCCORE
        7	Genome_NUCCORE
        8	Length_NUCCORE
        9	Source_NUCCORE
        10	UID_ASSEMBLY
        11	Status_ASSEMBLY
        12	SeqReleaseDate_ASSEMBLY
        13	SubmissionDate_ASSEMBLY
        14	Latest_ASSEMBLY
        15	UID_BIOSAMPLE
        16	ACC_BIOSAMPLE
        17	Location_BIOSAMPLE
        18	Coordinates_BIOSAMPLE
        19	IsolationSource_BIOSAMPLE
        20	Host_BIOSAMPLE
        21	SamplType_BIOSAMPLE
        22	taxon_name
        23	taxon_rank
        24	lineage
        25	taxon_species_id
        26	taxon_species_name
        27	taxon_genus_id
        28	taxon_genus_name
        29	taxon_family_id
        30	taxon_family_name
        31	taxon_order_id
        32	taxon_order_name
        33	taxon_class_id
        34	taxon_class_name
        35	taxon_phylum_id
        36	taxon_phylum_name
        37	taxon_superkingdom_id
        38	taxon_superkingdom_name
        39	loc_lat
        40	loc_lng
        41	loc_parsed
        42	GC_NUCCORE
        43	Identical
        44	OldVersion
        45	hits_rMLST
        46	hitscount_rMLST
        47	D1
        48	D2
        49	plasmidfinder
        50	pmlst
        '''
        with open(self.plsdb_info_table, 'r') as f:
            for n, row in enumerate(f):
                if n == 0:
                    continue 
                else:
                    data = row.rstrip().split('\t')
                    self.plsdb_entry2length[data[1]] = int(data[8])  



    def parse_plsdb(self,):
        
        # NODE_4	NZ_CP010410.1	91.139	79	2	5	17617	17693	454313	454238	1.27e-17	102	0 

        contig_regex = re.compile("(NODE_[0-9]+)_.*")

        for plsdb in self.plsdb:
            sample = plsdb.split("/")[1]
            self.sample2plsdb[sample] = {}   
            with open(plsdb) as f:
                for n, row in enumerate(f):
                    data = row.rstrip().split("\t")
                    contig_name = data[0]
                    plasmid_accession = data[1]
                    pident = data[2]
                    length = data[3]
                    mismatch = data[4]
                    gapopen = data[5]
                    qstart = data[6]
                    qend = data[7]
                    sstart = data[8]
                    send = data[9]
                    evalue = data[10]
                    bitscore = data[11]

                    plasmid_length = self.plsdb_entry2length[plasmid_accession] 
                    contig_length = int(self.sample2contig[sample][contig_name]["contig_length"])
                    plasmid_coverage = round((float(length) / plasmid_length) * 100, 2)
                    contig_coverage = round((float(length) / contig_length) * 100, 2)

                    if contig_coverage >= 80:
                        self.sample2contig[sample][contig_name]["plasmid"] = "TRUE" 
                    else:
                        self.sample2contig[sample][contig_name]["plasmid"] = "FALSE"

                    if contig_name not in self.sample2plsdb[sample]:
                        self.sample2plsdb[sample][contig_name] = [data + [plasmid_length, contig_length, plasmid_coverage, contig_coverage]] 




c = Contig(rgi, 
           core_genes, 
           plsdb, 
           COG_mobilome, 
           deepvirfinder, 
           depth, 
           GC,
           n_ORFs,
           plsdb_info_table)

c.parse_depth()
c.parse_GC()
c.parse_n_ORFs()
c.parse_deepvirfinder()
c.parse_cog_mobilome()
c.parse_core_genes()
c.parse_plsdb()
c.parse_rgi()

o = open(snakemake.output[1], 'w')

header =[
    "sample",
    "contig",
    "unique_id",
    "contig_length",
    "average_depth",
    "var_depth",
    "GC",
    "n_ORFs",
    "deepvirfinder_score",
    "deepvirfinder_pvalue",
    "plasmid",
    "COG_mobilome",
    "core_genes",  
    "rgi" 
] 


sql = 'create table contig (sample varchar(200), contig_name varchar(200),unique_id varchar(200), contig_length INTEGER, average_depth FLOAT, ' \
      ' var_depth FLOAT, GC FLOAT, n_ORFs INT, deepvirfinder_score FLOAT, deepvirfinder_pvalue FLOAT, plasmid boolean, ' \
      ' COG_mobilome INTEGER, CORE_genes INTEGER, rgi INTEGER, group_1 varchar(200), group_2 varchar(200))'

cursor.execute(sql,)

sql_template = 'insert into contig values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'

o.write('\t'.join(header) + '\n')
for sample in c.sample2contig:
    for contig in c.sample2contig[sample]: 

        unique_name = f"{sample}_{contig}"

        # {'contig_length': '500', 'average_depth': '2.54571', 'var_depth': '1.93343', 'GC': '65.2', 'deepvirfinder': {'score': '0.0006834288360551', 'pvalue': '0.8429979867918493'}, 'plasmid': 'FALSE'}
        contig_length =  c.sample2contig[sample][contig]["contig_length"] 
        average_depth =  c.sample2contig[sample][contig]["average_depth"] 
        var_depth =  c.sample2contig[sample][contig]["var_depth"] 
        GC =  c.sample2contig[sample][contig]["GC"] 
        try:
            n_ORFs =  c.sample2contig[sample][contig]["n_ORFs"] 
        except KeyError:
            n_ORFs = 0
        deepvirfinder_score =  c.sample2contig[sample][contig]["deepvirfinder"]["score"] 
        deepvirfinder_pvalue =  c.sample2contig[sample][contig]["deepvirfinder"]["pvalue"] 
        try:
            plasmid =  c.sample2contig[sample][contig]["plasmid"]
        except :
            plasmid = "FALSE"
        if plasmid == 'TRUE':
            plasmid = 1
        else:
            plasmid = 0
        try:
            COG_mobilome =  c.sample2contig[sample][contig]["COG_mobilome"]
        except :
            COG_mobilome = 0           
        try:
            core_genes =  c.sample2contig[sample][contig]["core_genes"]
        except :
            core_genes = 0 
        try:
            rgi =  c.sample2contig[sample][contig]["rgi"]
        except :
            rgi = 0
        
        group_1 = all_samples.loc[sample, "group_1"]
        group_2 = all_samples.loc[sample, "group_2"]

        lst = [sample, 
               contig,
               unique_name,
               contig_length,
               average_depth,
               var_depth,
               GC,
               n_ORFs,
               deepvirfinder_score,
               deepvirfinder_pvalue,
               plasmid,
               COG_mobilome,
               core_genes,
               rgi,
               group_1,
               group_2]
        cursor.execute(sql_template, lst)
        o.write('\t'.join([str(i) for i in (lst)]) + '\n') 

conn.commit()