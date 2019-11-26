
rule index_diamond_database:
    conda:
        "../../envs/diamond.yml"
    singularity:
        "docker://quay.io/biocontainers/diamond:0.9.28--h56fc30b_0"
    input:
        "reference_databases/{data_type}/{db_name}.faa"
    output:
        "reference_databases/{data_type}/{db_name}.dmnd"
    shell:
        """
        diamond makedb --in {input[0]} -d {output[0]}
        """

rule execute_diamond:
    conda:
        "../../envs/diamond.yml"
    singularity:
        "docker://quay.io/biocontainers/diamond:0.9.28--h56fc30b_0"
    input:
        "reference_databases/{data_type}/{db_name}.dmnd",
        "samples/{sample}/reads/raw/merged_combined.fastq"
    threads:
        8
    params:
        evalue = 0.01
    log:
        "samples/{sample}/{data_type}/{db_name}/diamond_merged_combined.log"
    output:
        "samples/{sample}/{data_type}/{db_name}/diamond_merged_combined.m8"
    shell:
        """
        db_path="{input[0]}"
        echo $db_path
        db_name=${{db_path%.*}}
        echo $db_name
        diamond blastx --evalue {params[0]} -p {threads} -d $db_name -q {input[1]} -o {output[0]} --max-target-seqs 1 --max-hsps 1 >> {log}
        """


