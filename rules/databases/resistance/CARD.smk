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
    shell:
        """
        tar -xvf {input[0]} -C reference_databases/resistance
        ar -xvzf {input[1]} -C reference_databases/resistance
        #mv reference_databases/resistance/card-data/* reference_databases/resistance/
        mv reference_databases/resistance/protein_fasta_protein_homolog_model.fasta reference_databases/resistance/CARD_protein_homolog_model.faa
        """
