#!/usr/bin/env python

input_fasta_list = snakemake.input.fasta_files
output_file = snakemake.output[0]


from Bio import SeqIO
import sys
import os
import re

o = open(output_file, 'w')
for fasta in input_fasta_list:
    with open(fasta, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            o.write("%s\t%s\n" % (record.name,
                                  re.sub('\.','_',os.path.basename(fasta.split(".fasta")[0]))))
