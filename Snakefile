
include: "rules/making_sample_dataset.rules"

print()

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


rule homology_search:
   input:
       expand("samples/{sample}/mmseq_search/virulence/VF_db/best_hits.m8", sample = list(read_naming.keys())),
       "report/mmseq_search/virulence/VF_db/RPKM.db"

pipeline_path = workflow.basedir + '/'
multiqc_configfile = pipeline_path + "data/configuration_files/multiqc/config.yaml"



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
include: "rules/QC/decontam_reads.rules"

# metagenomic assembly
include: "rules/assembly_and_binning/quality_check.rules"
include: "rules/assembly_and_binning/prodigal.rules"
include: "rules/assembly_and_binning/bwa.rules"
include: "rules/assembly_and_binning/assembly.rules"
include: "rules/assembly_and_binning/binning.rules"

# database setup
include: "rules/databases/virulence/virulence.smk"
include: "rules/databases/mmseqs2.smk"

# taxonomy abundance
include: "rules/taxonomy_abundance/motus2.rules"
include: "rules/protein_abundance/calculate_RPKM.rules"

#taxonomy profiling
#kaiju
include: "rules/taxonomy_profiling/kaiju.rules"
include: "rules/taxonomy_profiling/kraken2.rules"