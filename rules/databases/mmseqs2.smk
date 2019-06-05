
rule index_mmseqs2_database:
    conda:
        "../../envs/mmseqs.smk"
    input:
        "reference_databases/{data_type}/{db_name}.faa"
    output:
        "reference_databases/{data_type}/{db_name}_mmseqsDB"
    shell:
        """
        mmseqs createdb {input[0]} {output[0]}
        """


rule execute_mmseqs2:
    conda:
        "../../envs/mmseqs.smk"
    input:
        "reference_databases/{data_type}/{db_name}_mmseqsDB",
        "samples/{sample}/reads/trimmed/{pair}_paired.fastq"
    threads:
        8
    params:
        min_seq_id = 0.9
    output:
        "samples/{sample}/mmseq_search/{data_type}/{db_name}/{pair}.m8"
    shell:
        """
        mmseqs easy-search --threads {threads} --min-seq-id {params[0]} {input[1]} {input[0]} {output[0]} tmp
        """

rule extract_best_hit:
    shell:
       """
       sort -u -k1,1 --merge samples/B6T_subset/mmseq_search/virulence/VF_db/R1.m8
       """
