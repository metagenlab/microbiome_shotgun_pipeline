

rule prepare_virulence_fasta:
    conda:
        "../../../envs/biopython.yml"
    input:
        config["virulence_database"],
        config["virulence_fasta"]
    output:
        "reference_databases/virulence/VF_db.faa"
    script: "scripts/extract_database_from_sqlitedb.py"
