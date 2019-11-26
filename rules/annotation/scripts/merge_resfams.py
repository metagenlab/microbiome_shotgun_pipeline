
import pandas

resfams_files = snakemake.input
sample_table = snakemake.params[0]

all_samples = pandas.read_csv(sample_table, sep="\t", index_col=0)

def parse_hmmsearch(hmmsearch_result):
    from Bio import SearchIO

    hmmer_records = [i for i in SearchIO.parse(hmmsearch_result, 'hmmer3-text')]

    hsp_list = []

    for hmmer_record in hmmer_records:

        try:
            best_hit_id = hmmer_record.hits[0].id
        except IndexError:
            continue
        
        profile_id =  hmmer_record.id

        '''
        'append', 'bias', 'bitscore', 'description', 'description_all', 'domain_exp_num', 'domain_obs_num',
        'evalue', 'filter', 'fragments', 'hsps', 'id',
        'id_all', 'index', 'is_included', 'map', 'pop', 'query_description', 'query_id', 'sort'

        '''
        for hit in hmmer_record.hits:
            for hsp in hit.hsps:
                print(hsp)
                result = {
                    "profile_id" : hmmer_record.accession,
                    "profile_length" : hmmer_record.seq_len,
                    "best_hit_id" : hit.id,
                    "bias" : hit.bias,
                    "bitscore" : hit.bitscore,
                    "hit_start" : hsp.hit_start,
                    "hit_end" : hsp.hit_end,
                    "query_start" : hsp.query_start,
                    "query_end" : hsp.query_end,
                    "evalue" : hit.evalue,
                    "query_coverage" : round((hsp.query_end-hsp.query_start)/float(hmmer_record.seq_len),2)
                    }
                hsp_list.append(result)

    return hsp_list

header1 = ["tool",
           "sample", 
           "group_1",
           "group_2", 
           "contig_id",
           "ORF_id",
           "profile_id",
           "profile_length",
           "bias",
           "bitscore",
           "hit_start",
           "hit_end",
           "query_start",
           "query_end",
           "evalue",
           "query_coverage"]

o_resfams_hits = open(snakemake.output[0], "w")
o_resfams_hits.write("\t".join(header1) + '\n')

for resfams_file in resfams_files:
    hsp_list = parse_hmmsearch(resfams_file)
    sample_name = resfams_file.split("/")[1]
    group_1 = all_samples.loc[sample_name, "group_1"]
    group_2 = all_samples.loc[sample_name, "group_2"]


    for hsp in hsp_list:
        if "plass" in resfams_file:
            tool = "plass"
            contig_id = None
            orf_number = hsp["best_hit_id"] #.split("[Orf: ")[1].split(",")[0]
        if "prodigal" in resfams_file:
            tool = "prodigal"
            contig_prefix, contig_number, orf_number = hsp["best_hit_id"].split(" ")[0].split("_")
            contig_id = "%s_%s" % (contig_prefix, contig_number)
        o_resfams_hits.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tool,
                                                                                                   sample_name,
                                                                                                   group_1,
                                                                                                   group_2,
                                                                                                   contig_id,
                                                                                                   orf_number,
                                                                                                   hsp["profile_id"],
                                                                                                   hsp["profile_length"],
                                                                                                   hsp["bias"],
                                                                                                   hsp["bitscore"],
                                                                                                   hsp["hit_start"],
                                                                                                   hsp["hit_end"],
                                                                                                   hsp["query_start"],
                                                                                                   hsp["query_end"],
                                                                                                   hsp["evalue"],
                                                                                                   hsp["query_coverage"]
                                                                                                   ))