
rule index_mmseqs2_database:
    conda:
        "../../envs/mmseqs.yml"
    singularity:
        "docker://quay.io/biocontainers/mmseqs2:10.6d92c--h2d02072_0"
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
        "../../envs/mmseqs.yml"
    singularity:
        "docker://quay.io/biocontainers/mmseqs2:10.6d92c--h2d02072_0"
    input:
        "reference_databases/{data_type}/{db_name}_mmseqsDB",
        "samples/{sample}/reads/trimmed/{pair}_paired.fastq"
    threads:
        8
    params:
        min_seq_id = 0.9,
        min_aln_len = 24
    output:
        "samples/{sample}/mmseq_search/{data_type}/{db_name}/{pair}.m8"
    shell:
        """
        mmseqs easy-search --threads {threads} --min-seq-id {params[0]} --min-aln-len {params[1]} {input[1]} {input[0]} {output[0]} tmp
        """

rule extract_best_hit:
    input:
        "samples/{sample}/mmseq_search/{data_type}/{db_name}/R1.m8",
        "samples/{sample}/mmseq_search/{data_type}/{db_name}/R2.m8"
    output:
        "samples/{sample}/mmseq_search/{data_type}/{db_name}/best_hits.m8"
    shell:
       """
       sort -u -k1,1 --merge {input[0]} > {output[0]}
       sort -u -k1,1 --merge {input[1]} >> {output[0]}
       """
