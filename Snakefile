
include:
    "rules/making_sample_dataset.rules"

rule all:
    input:
       qualimap_report = "report/multiqc_assembly/multiqc_report.html",
       buble_plot_maxbin = expand("samples/{sample}/binning/maxbin/gc_cov_buble_test.svg", sample = list(read_naming.keys())),
       kaiju_phylum = expand("samples/{sample}/binning/kaiju/report_phylum.txt", sample = list(read_naming.keys())),
       #chekm_maxbin = expand("samples/{sample}/binning/maxbin/checkm_bacteria/storage/bin_stats.analyze.tsv", sample = list(read_naming.keys()))


pipeline_path = workflow.basedir + '/'
multiqc_configfile = pipeline_path + "data/configuration_files/multiqc/config.yaml"

include:
    "rules/assembly.rules"
include:
    "rules/binning.rules"
include:
    "rules/logging.rules"
include:
    "rules/quality_check.rules"
include:
    "rules/trimmomatic.rules"
include:
    "rules/prodigal.rules"
include:
    "rules/get_reads.rules"
include:
    "rules/get_sras.rules"
include:
    "rules/multiqc.rules"
include:
    "rules/bwa.rules"
include:
    "rules/qualimap.rules"
include:
    "rules/quast.rules"
include:
    "rules/fastqc.rules"
