
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
        "samples/{sample}/reads/raw/merged_combined.fastq"
    threads:
        8
    params:
        evalue = 0.001,
    log:
        "samples/{sample}/{data_type}/{db_name}/mmseqs_merged_combined.log"
    output:
        "samples/{sample}/{data_type}/{db_name}/mmseqs_merged_combined_raw.m8"
    shell:
        """
        mmseqs easy-search --threads {threads} -e {params[0]} {input[1]} {input[0]} {output[0]} $(dirname {output[0]})/tmp >> {log}
        """


rule update_identity_values:
    input:
        "samples/{sample}/{data_type}/{db_name}/mmseqs_merged_combined_raw.m8"
    threads:
        1
    output:
        "samples/{sample}/{data_type}/{db_name}/mmseqs_merged_combined.m8"
    shell:
        """
        awk '{{$3=sprintf("%.1f",$3*100)}}1' OFS="\t" {input[0]} >> {output[0]}
        """