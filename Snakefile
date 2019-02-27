
rule all:
    input:
       qualimap_report = "quality/multiqc/assembly/multiqc_report.html",
       buble_plot_maxbin = expand("samples/{sample}/binning/maxbin/gc_cov_buble_test.svg", sample = list(read_naming.keys())),
       kaiju_phylum = expand("samples/{sample}/binning/kaiju/report_phylum.txt", sample = list(read_naming.keys())),
       #chekm_maxbin = expand("samples/{sample}/binning/maxbin/checkm_bacteria/storage/bin_stats.analyze.tsv", sample = list(read_naming.keys()))

inlude:
    "rules/assembly.rules"
inlude:
    "rules/binning.rules"
inlude:
    "rules/logging.rules"
inlude:
    "rules/making_sample_dataset.rules"
inlude:
    "rules/quality_check.rules"
inlude:
    "rules/trimmomatic.rules"
