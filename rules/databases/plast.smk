
rule execute_plast:
    conda:
        "../../envs/plast.yml"
    singularity:
        "docker://metagenlab/plast:2.3.2"
    input:
        "reference_databases/{data_type}/{db_name}.faa",
        "samples/{sample}/reads/raw/merged_combined.fasta"
    threads:
        8
    params:
        evalue = 0.01
    log:
        "samples/{sample}/{data_type}/{db_name}/plast_merged_combined.log"
    output:
        "samples/{sample}/{data_type}/{db_name}/plast_merged_combined.m8"
    shell:
        """
        plast -p plastx -a {threads} -d {input[0]} -i {input[1]} -M BLOSUM62 -s 45 -seeds-use-ratio 60 -max-database-size 50000000 -e {params[0]} -G 11 -E 1 -o {output[0]} -F F -bargraph -verbose -force-query-order 1000 -max-hit-per-query 100 -max-hsp-per-hit 1 >> {log};
        """
