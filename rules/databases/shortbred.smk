


rule execute_shortbred:
    conda:
        "../../envs/plast.yml"
    singularity:
        "docker://metagenlab/shortbred:1.0"
    input:
        "reference_databases/{data_type}/shortbred/{db_name}.faa",
        "samples/{sample}/reads/raw/merged_combined.fastq"
    threads:
        12
    params:
        evalue = 0.01
    log:
        "samples/{sample}/{data_type}/{db_name}/shortbred_merged_combined.log"
    output:
        "samples/{sample}/{data_type}/{db_name}/shortbred_merged_combined.txt"
    shell:
        """
        echo shortbred_quantify.py --markers {input[0]}  --wgs {input[1]}  --results {output[0]}  --tmp  $(dirname {output[0]})/sortbred_quantify --threads {threads} 
        shortbred_quantify.py --markers {input[0]}  --wgs {input[1]}  --results {output[0]}  --tmp  $(dirname {output[0]})/sortbred_quantify --threads {threads}
        """

rule all_shortbred:
    conda:
        "../../envs/python-r.yml"
    singularity:
        "docker://metagenlab/microbiome-shotgun-pipeline:1.0"
    params:
        sample_table = config["local_samples"]
    input:
        expand("samples/{sample}/{{data_type}}/{{db_name}}/shortbred_merged_combined.txt", sample=list(read_naming.keys()))
    output:
        "report/{data_type}/{db_name}/shortbred.tab",
    script:
        "scripts/combine_shortbred.py"