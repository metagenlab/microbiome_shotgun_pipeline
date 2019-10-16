

rule prepare_virulence_fasta:
    conda:
        "../../envs/python-r.yml"
    singularity:
        "docker://metagenlab/microbiome-shotgun-pipeline:1.0"
    input:
        config["virulence_database"],
        config["virulence_fasta"]
    output:
        "reference_databases/virulence/VF_db.faa"
    script: "scripts/extract_database_from_sqlitedb.py"
