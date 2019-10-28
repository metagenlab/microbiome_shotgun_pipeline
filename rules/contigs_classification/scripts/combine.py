#!/usr/bin/env python


rgi = snakemake.input["rgi"]
core_genes = snakemake.input["core_genes"]
plsdb = snakemake.input["plsdb"]
COG_mobilome = snakemake.input["COG_mobilome"]
deepvirfinder = snakemake.input["deepvirfinder"]
depth = snakemake.input["depth"]
cat_taxonomy = snakemake.input["cat_taxonomy"]
GC = snakemake.input["GC"]

from Bio import SearchIO

class Contig:

    def __init__(rgi, 
                 core_genes, 
                 plsdb, 
                 COG_mobilome, 
                 deepvirfinder, 
                 depth, 
                 cat_taxonomy, 
                 GC):

        self.rgi = rgi 
        self.core_genes = core_genes # ok
        self.plsdb = plsdb
        self.COG_mobilome = COG_mobilome # ok
        self.deepvirfinder = deepvirfinder # ok
        self.depth = depth # ok
        self.cat_taxonomy = cat_taxonomy
        self.GC = GC # ok

        self.sample2contig = {}

def parse_depth(self):

    '''
    return dictionnary contig 2 len, average depth, var depth
    '''
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
                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    self.sample2contig[sample][contig_name]["GC"] = gc


    def parse_deepvirfinder():
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
    

    def parse_cog_mobilome():
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

        contig_regex = re.compile("(NODE_[0-9]+)_.*")

        for COG_mobilome in self.COG_mobilome:
            sample = COG_mobilome.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
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

                    contig_name = contig_regex.search(ORF_id).group(0)

                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    if "COG_mobilome" not in self.sample2contig[sample][contig_name]:
                        self.sample2contig[sample][contig_name]["COG_mobilome"] = {}
                    # only keep first hit
                    if ORF_id not in self.sample2contig[sample][contig_name]["COG_mobilome"]:
                        self.sample2contig[sample][contig_name]["COG_mobilome"][ORF_id] = data[1:len(data)]  


    def parse_core_genes():
        contig_regex = re.compile("(NODE_[0-9]+)_.*")

        for core_genes in self.core_genes:
            sample = core_genes.split("/")[1]
            if sample not in self.sample2contig:
                self.sample2contig[sample] = {}
            records = SearchIO.parse(core_genes, 'hmmer3-tab')
            for record in records:
                for hit in record.hits:
                    ORF_id = hit.id
                    HMM_id = record.id
                    evalue = hit.evalue
                    contig_name = contig_regex.search(ORF_id).group(0)

                    if contig_name not in self.sample2contig[sample]:
                        self.sample2contig[sample][contig_name] = {}
                    if "core_genes" not in self.sample2contig[sample][contig_name]:
                        self.sample2contig[sample][contig_name]["core_genes"] = {}
                    self.sample2contig[sample][contig_name]["core_genes"][ORF_id] = [HMM_id, evalue] 

    def parse_rgi():

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

