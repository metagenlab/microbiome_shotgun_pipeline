
rule calculate_depth_from_reads_hits:
    conda:
        "../../envs/python-r.yml"
    singularity:
        "docker://metagenlab/diag-pipeline-python-r:1.1"
    input:
        "samples/{sample}/{data_type}/{db_name}/{search_tool}_best_hits_{identity}.m8",
        "reference_databases/{data_type}/{db_name}.faa"
    output:
        "samples/{sample}/{data_type}/{db_name}/{search_tool}_best_hits_{identity}_depth_plots.pdf",
        "samples/{sample}/{data_type}/{db_name}/{search_tool}_best_hits_{identity}_statistics_per_gene.tab",
        "samples/{sample}/{data_type}/{db_name}/{search_tool}_best_hits_{identity}_depth_per_gene.tab",
    script: "scripts/calculate_depth_from_tabulated_blast.py"


rule all_gene_statistics:
    conda:
        "../../envs/python-r.yml"
    singularity:
        "docker://metagenlab/microbiome-shotgun-pipeline:1.0"
    params:
        sample_table = config["local_samples"],
    input:
        statistic_files = expand("samples/{sample}/{{data_type}}/{{db_name}}/{{search_tool}}_best_hits_{{identity}}_statistics_per_gene.tab", sample=list(read_naming.keys())),
        silix_98 = "reference_databases/{data_type}/{db_name}_0.98_0.6.fnodes",
        silix_95 = "reference_databases/{data_type}/{db_name}_0.95_0.6.fnodes",
        silix_90 = "reference_databases/{data_type}/{db_name}_0.90_0.6.fnodes",
        silix_80 = "reference_databases/{data_type}/{db_name}_0.80_0.6.fnodes",
    output:
        "report/{data_type}/{db_name}/{search_tool}_best_hits_{identity}_statistics_per_gene.tab",
    script:
        "scripts/combine_gene_statistics.py"