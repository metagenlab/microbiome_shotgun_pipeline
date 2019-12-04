rule download_CARD_ontology:
    output:
        "reference_databases/resistance/ontology-v3.0.3.tar.gz",
    log:
        "reference_databases/resistance/download.log",
    shell:
        """
        curl https://card.mcmaster.ca/download/5/ontology-v3.0.3.tar.gz -o {output[0]} >> {log}
        """

rule download_CARD_data:
    output:
        "reference_databases/resistance/broadstreet-v3.0.3.tar.gz",
    log:
        "reference_databases/resistance/download.log",
    shell:
        """
        curl https://card.mcmaster.ca/download/0/broadstreet-v3.0.3.tar.gz -o {output[0]} >> {log}
        """

rule extract_CARD:
    input:
        "reference_databases/resistance/ontology-v3.0.3.tar.gz",
        "reference_databases/resistance/broadstreet-v3.0.3.tar.gz",
    output:
        "reference_databases/resistance/CARD_protein_homolog_model.faa",
        "reference_databases/resistance/card.json",
        "reference_databases/resistance/aro_index.tsv",
        "reference_databases/resistance/aro.tsv",
        "reference_databases/resistance/aro_categories_index.tsv",
    shell:
        """
        tar -xvf {input[0]} -C reference_databases/resistance
        tar -xvf {input[1]} -C reference_databases/resistance
        mv reference_databases/resistance/protein_fasta_protein_homolog_model.fasta reference_databases/resistance/CARD_protein_homolog_model.faa
        """

rule format_aro:
    input:
        "reference_databases/resistance/aro_index.tsv",
        "reference_databases/resistance/aro_categories_index.tsv",
    output:
        "reference_databases/resistance/accession2AMR_family.tab",
        "reference_databases/resistance/accession2drug_class.tab",
        "reference_databases/resistance/accession2resistance_mechanism.tab",
    script:
        "scripts/aro2data.py"