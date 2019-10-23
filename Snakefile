
include: "rules/making_sample_dataset.rules"

rule assembly:
    input:
       qualimap_report = "report/multiqc_assembly/multiqc_report.html",
       buble_plot_maxbin = expand("samples/{sample}/binning/maxbin/gc_cov_buble_test.svg", sample = list(read_naming.keys())),
       kaiju_phylum = expand("samples/{sample}/binning/kaiju/report_phylum.txt", sample = list(read_naming.keys())),
       #chekm_maxbin = expand("samples/{sample}/binning/maxbin/checkm_bacteria/storage/bin_stats.analyze.tsv", sample = list(read_naming.keys()))

rule databases_setup:
   input:
       "reference_databases/virulence/VF_db.faa",
       "reference_databases/virulence/VF_db_mmseqsDB"

rule motus:
    input:
        "report/motus2/relative/phylum/merged_profile.motus",
        "report/motus2/relative/genus/merged_profile.motus",
        "report/motus2/relative/species/merged_profile.motus",
        "report/motus2/relative/mOTU/merged_profile.motus",

rule homology_search:
   input:
       expand("samples/{sample}/mmseq_search/virulence/VF_db/best_hits.m8", sample = list(read_naming.keys())),
       "report/mmseq_search/virulence/VF_db/annotation.log",
       expand("samples/{sample}/mmseq_search/resistance/CARD_protein_homolog_model/best_hits.m8", sample = list(read_naming.keys())),
       "report/mmseq_search/resistance/CARD_protein_homolog_model/annotation.log",
       "report/mmseq_search/virulence/VF_db/plots/VF_counts.svg",
       "report/mmseq_search/resistance/CARD_protein_homolog_model/plots/RES_counts.svg"


rule anvio:
    input:
        expand("samples/{sample}/assembly/spades/large_contigs_edit.fasta", sample = list(read_naming.keys())),
        expand("samples/{sample}/gene_call/{sample}.faa", sample = list(read_naming.keys())),
        #expand("samples/{sample}/anvio/contigs_db/contigs.db", sample = list(read_naming.keys())),

rule annotation:
    input:
        expand("samples/{sample}/annotation/rgi/rgi.json", sample = list(read_naming.keys())),
        expand("samples/{sample}/contigs_classification/deepvirfinder/large_contigs_edit.fasta_gt500bp_dvfpred.txt", sample = list(read_naming.keys()))


pipeline_path = workflow.basedir + '/'
multiqc_configfile = "../../data/configuration_files/multiqc/config.yaml"

# general
include: "rules/logging.rules"

# read manipulation
include: "rules/read_manipulation/get_reads.rules"
include: "rules/read_manipulation/get_sras.rules"

# quality check
include: "rules/QC/multiqc.rules"
include: "rules/QC/qualimap.rules"
include: "rules/QC/quast.rules"
include: "rules/QC/fastqc.rules"
include: "rules/QC/trimmomatic.rules"

# metagenomic assembly
include: "rules/assembly_and_binning/quality_check.rules"
include: "rules/assembly_and_binning/prodigal.rules"
include: "rules/assembly_and_binning/bwa.rules"
include: "rules/assembly_and_binning/assembly.rules"
include: "rules/assembly_and_binning/binning.rules"

# database setup
include: "rules/databases/virulence/virulence.smk"
include: "rules/databases/resistance/CARD.smk"
include: "rules/databases/mmseqs2.smk"

# taxonomy abundance
include: "rules/taxonomy_abundance/motus2.rules"
include: "rules/protein_abundance/calculate_RPKM.rules"
include: "rules/protein_abundance/visualization.rules"

# anvio 
include: "rules/anvio/anvio.rules"

# annotation
include: "rules/annotation/rgi.rules"

# contig classification
include: "rules/contigs_classification/deepvirfinder.rules"
