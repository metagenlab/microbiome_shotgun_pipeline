

from Bio import SeqIO

def db2gene_size(db_fasta, split_id=False):
    gene2size = {}
    records = SeqIO.parse(db_fasta, "fasta")
    for record in records:
        if split_id:
            data = record.id.split("|")
            if len(data) == 4:
                name = data[1]
            else:
                name = record.id
        else:
            name = record.id

        gene2size[name] = len(record.seq)
    return gene2size

def get_gene_depth(blast_file, db_file, split_id):
    gene2size = db2gene_size(db_file, split_id=split_id)
    gene2depth = {}
    gene2n_hits = {}
    with open(blast_file, 'r') as f:
        for row in f:
            data = row.split("\t")
            hit_name = data[1]
            # python index stats from 0
            hit_start = int(data[8]) - 1
            hit_end = int(data[9])
            if hit_name not in gene2depth:
                # initiallize empty count list
                try:
                    gene2depth[hit_name] = [0] * gene2size[hit_name]
                    gene2n_hits[hit_name] = 0
                except:
                    print("malformed fasta entry: %s -- skipping" % (hit_name))
                    continue

            # increment covered region
            gene2depth[hit_name][hit_start:hit_end] = [i+1 for i in gene2depth[hit_name][hit_start:hit_end]]
            gene2n_hits[hit_name] += 1
    return gene2depth, gene2n_hits

def get_depth_and_coverage(gene2depth, 
                           gene2n_hits,
                           pdf_plots = "multipage_pdf.pdf",
                           statistics_per_gene = "statistics_per_gene.tab",
                           depth_per_gene = "depth_per_gene.tab"):
    import statistics
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    pdf = PdfPages(pdf_plots)
    o = open(statistics_per_gene, "w")
    p = open(depth_per_gene, "w")

    for n, gene in enumerate(gene2depth):
        if n == 0:
            o.write("accession\tARO\tgene_name\tn_hits\tmedian_depth\tgene_length\tpercent_covered\n")
        # gb|ACJ41739.1|ARO:3000780|adeI
        if "ARO" in gene and "|" in gene:
            data = gene.split("|")
            accession = data[1]
            aro_accession = data[2]
            gene_name = data[3]
        else:
             accession = gene
             aro_accession = '-'
             gene_name = '-'

        n_covered = sum(1 for i in gene2depth[gene] if i > 0)
        gene_length = len(gene2depth[gene])
        percent_covered = round((n_covered/float(gene_length)) * 100, 2)
        median_depth = statistics.median(gene2depth[gene])
        o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (accession, 
                                                  aro_accession,
                                                  gene_name,
                                                  gene2n_hits[gene], 
                                                  median_depth,
                                                  gene_length,
                                                  percent_covered))
        if percent_covered > 30 and median_depth > 3:
            plt.plot(gene2depth[gene])
            plt.ylabel(gene)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

        p.write("%s\t%s\n" % (gene, gene2depth[gene]))
    pdf.close()

input_m8 = snakemake.input[0]
database_fasta = snakemake.input[1]
plots = snakemake.output[0]
statistics = snakemake.output[1]
details = snakemake.output[2]

if 'mmseqs' in input_m8:
    gene2depth, gene2n_hits = get_gene_depth(input_m8, database_fasta, split_id=True)
else:
    gene2depth, gene2n_hits = get_gene_depth(input_m8, database_fasta, False)

get_depth_and_coverage(gene2depth, 
                       gene2n_hits,
                       plots,
                       statistics,
                       details)
