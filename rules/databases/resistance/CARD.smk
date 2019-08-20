rule download_CARD:
    output:
        "reference_databases/resistance/ontology-v3.0.3.tar.gz",
    log:
        "reference_databases/resistance/download.log",
    shell:
        """
        curl https://card.mcmaster.ca/download/5/ontology-v3.0.3.tar.gz -o {output[0]} >> {log}
        """

rule extract_CARD:
    input:
        "reference_databases/resistance/ontology-v3.0.3.tar.gz",
    output:
        "reference_databases/resistance/CARD_protein_homolog_model.faa",
        "reference_databases/resistance/card.json",
    shell:
        """
        tar -xvzf {input[0]} -C reference_databases/resistance
        mv reference_databases/resistance/card-data/* reference_databases/resistance/
        mv reference_databases/resistance/protein_fasta_protein_homolog_model.fasta reference_databases/resistance/CARD_protein_homolog_model.faa
        """
