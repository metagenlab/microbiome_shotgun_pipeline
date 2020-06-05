import pandas as pd
from Bio import SeqIO
blastntab=snakemake.input[0]#tabular blastn output
min_identity=snakemake.params.ident
contigs=snakemake.params.contigs
out=snakemake.output[0]
try:
    tab=pd.read_csv(blastntab,sep='\t',header=None)
    tab.columns = ['query_id', 'subject_id ', 'pct_identity', 'aln_length', 'n_of_mismatches', 'gap_openings', 'q_start',
             'q_end', 's_start', 's_end', 'e_value', 'bit_score']

    if len(tab) == 1:
        contig_headers = list(tab['query_id'])
    else:
        contig_headers = list(set(tab[tab.pct_identity>=min_identity]['query_id']))

    handle = SeqIO.parse(contigs, 'fasta')
    records=[]
    for record in handle:
        if record.id in contig_headers:
            records.append(record)
    SeqIO.write(records, out, 'fasta')
except ValueError:
     f=open(out,'w')
     f.write('>no_contigs_found')
     f.close()