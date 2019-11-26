
rule extract_best_hit:
    input:
        "samples/{sample}/{data_type}/{db_name}/{search_tool}_merged_combined.m8",
    output:
        "samples/{sample}/{data_type}/{db_name}/{search_tool}_best_hits.m8"
    shell:
       """
       sort -u -k1,1 --merge {input[0]} > {output[0]}
       """


rule filter_based_on_identity:
    input:
        "samples/{sample}/{data_type}/{db_name}/{search_tool}_best_hits.m8",
    output:
        "samples/{sample}/{data_type}/{db_name}/{search_tool}_best_hits_{identity}.m8"
    shell:
       """
       cat {input[0]} | awk '$3>{wildcards.identity}' > {output[0]}
       """