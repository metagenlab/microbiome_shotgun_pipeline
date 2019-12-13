
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



def get_mummer_lst():
    import itertools
    lst = list(read_naming.keys())
    combinations = [i for i in itertools.combinations(lst, 2)]
    out = []
    for pair in combinations: 
        out.append("report/alignments/mummer/nucmer_%s_vs_%s.delta" % (pair[0], pair[1]))
    return out 


rule annotation:
    input:
        # expand("samples/{sample}/annotation/rgi/rgi_proteins.json", sample = list(read_naming.keys())),
        #expand("samples/{sample}/contigs_classification/core_genes/{sample}.tsv", sample = list(read_naming.keys())),
        #expand("samples/{sample}/contigs_classification/plsdb/{sample}.tsv", sample = list(read_naming.keys())),
        #expand("samples/{sample}/contigs_classification/COG_mobilome/{sample}.tsv", sample = list(read_naming.keys())),
        #expand("samples/{sample}/contigs_classification/deepvirfinder/large_contigs_edit.fasta_gt500bp_dvfpred.txt", sample = list(read_naming.keys())),
        #expand("samples/{sample}/bwa/{sample}_assembled/{sample}_metabat2_depth.txt", sample = list(read_naming.keys())),
        #expand("samples/{sample}/binning/kaiju/report_phylum.txt", sample = list(read_naming.keys())),
        # "samples/{sample}/contigs_classification/cat_bat/{tool}_ORFs_out.CAT.contig2classification.txt"
        expand("samples/{sample}/contigs_classification/cat_bat/prodigal_ORFs.CAT.contig2classification_names.txt", sample = list(read_naming.keys())),
        expand("samples/{sample}/contigs_classification/cat_bat/contigs.CAT.contig2classification_names.txt", sample = list(read_naming.keys())),
        expand("samples/{sample}/contigs_classification/cat_bat/prodigal_ORFs.CAT.summary", sample = list(read_naming.keys())),
        expand("samples/{sample}/contigs_classification/cat_bat/contigs.CAT.summary", sample = list(read_naming.keys())),


        "report/resistance/CARD_protein_homolog_model/diamond_RPKM_98_annotated.tsv",
        "report/resistance/CARD_protein_homolog_model_edit/mmseqs_RPKM_98_annotated.tsv",
        "report/resistance/CARD_protein_homolog_model_edit/plast_RPKM_98_annotated.tsv",
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit/diamond_best_hits_98_depth_plots.pdf", sample = list(read_naming.keys())),      
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit/plast_best_hits_98_depth_plots.pdf", sample = list(read_naming.keys())),
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit/plast_best_hits_98_depth_plots.pdf", sample = list(read_naming.keys())),
        "report/resistance/CARD_protein_homolog_model/diamond_best_hits_98_statistics_per_gene.tab",
        "report/resistance/CARD_protein_homolog_model/mmseqs_best_hits_98_statistics_per_gene.tab",
        "report/resistance/CARD_protein_homolog_model/plast_best_hits_98_statistics_per_gene.tab",
        

        "report/resistance/CARD_protein_homolog_model_edit/diamond_RPKM_90_annotated.tsv",
        "report/resistance/CARD_protein_homolog_model_edit/mmseqs_RPKM_90_annotated.tsv",
        "report/resistance/CARD_protein_homolog_model_edit/plast_RPKM_90_annotated.tsv",
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit/diamond_best_hits_90_depth_plots.pdf", sample = list(read_naming.keys())),
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit//mmseqs_best_hits_90_depth_plots.pdf", sample = list(read_naming.keys())),
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit/plast_best_hits_90_depth_plots.pdf", sample = list(read_naming.keys())),
        "report/resistance/CARD_protein_homolog_model/diamond_best_hits_90_statistics_per_gene.tab",
        "report/resistance/CARD_protein_homolog_model/mmseqs_best_hits_90_statistics_per_gene.tab",
        "report/resistance/CARD_protein_homolog_model/plast_best_hits_90_statistics_per_gene.tab",

        "report/resistance/CARD_protein_homolog_model_edit/diamond_RPKM_80_annotated.tsv",
        "report/resistance/CARD_protein_homolog_model_edit/mmseqs_RPKM_80_annotated.tsv",
        "report/resistance/CARD_protein_homolog_model_edit/plast_RPKM_80_annotated.tsv",
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit/diamond_best_hits_80_depth_plots.pdf", sample = list(read_naming.keys())),
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit/plast_best_hits_80_depth_plots.pdf", sample = list(read_naming.keys())),
        expand("samples/{sample}/resistance/CARD_protein_homolog_model_edit/plast_best_hits_80_depth_plots.pdf", sample = list(read_naming.keys())),
        "report/resistance/CARD_protein_homolog_model/diamond_best_hits_80_statistics_per_gene.tab",
        "report/resistance/CARD_protein_homolog_model/mmseqs_best_hits_80_statistics_per_gene.tab",
        "report/resistance/CARD_protein_homolog_model/plast_best_hits_80_statistics_per_gene.tab",


        "report/annotations/rgi/prodigal_ORFs_rgi.tab",
        "report/annotations/rgi/plass_ORFs_rgi.tab",
        "report/contigs_classification/summary.tsv",
        "reference_databases/resistance/accession2AMR_family.tab",
        "report/annotations/rgi/contigs_rgi.tab",

        # "report/resistance/CARD_shortbred_markers/shortbred.tab",
        # "report/annotations/rgi/prodigal_ORFs_amrfinderplus.tab",
        # "report/annotations/rgi/plass_ORFs_amrfinderplus.tab",
        get_mummer_lst(),
        #expand("samples/{sample}/contigs_classification/{sample}_GC.tsv", sample = list(read_naming.keys())),


pipeline_path = workflow.basedir + '/'
multiqc_configfile = "../../data/configuration_files/multiqc/config.yaml"

# general
include: "rules/logging.rules"

# db clustering 

include: "rules/clustering/silix.rules"

# read manipulation
include: "rules/read_manipulation/get_reads.rules"
include: "rules/read_manipulation/get_sras.rules"
include: "rules/read_manipulation/merge_reads.rules"

# quality check
include: "rules/QC/multiqc.rules"
include: "rules/QC/qualimap.rules"
include: "rules/QC/quast.rules"
include: "rules/QC/fastqc.rules"
include: "rules/QC/trimmomatic.rules"

# metagenomic assembly
include: "rules/assembly_and_binning/gc_cov_plots.rules"
include: "rules/assembly_and_binning/prodigal.rules"
include: "rules/assembly_and_binning/plass.rules"
include: "rules/assembly_and_binning/bwa.rules"
include: "rules/assembly_and_binning/assembly.rules"
include: "rules/assembly_and_binning/binning.rules"
include: "rules/assembly_and_binning/checkm.rules"
include: "rules/assembly_and_binning/statistics.rules"

# database setup
include: "rules/databases/virulence/virulence.smk"
include: "rules/databases/resistance/CARD.smk"
include: "rules/databases/mmseqs2.smk"
include: "rules/databases/diamond.smk"
include: "rules/databases/plast.smk"
include: "rules/databases/filter_hits.smk"
include: "rules/databases/genes_statistics.smk"
include: "rules/databases/shortbred.smk"

# taxonomy abundance
include: "rules/taxonomy_abundance/motus2.rules"
include: "rules/protein_abundance/calculate_RPKM.rules"
include: "rules/protein_abundance/visualization.rules"

# anvio 
include: "rules/anvio/anvio.rules"

# annotation
include: "rules/annotation/rgi.rules"
include: "rules/annotation/amrfinder.rules"
include: "rules/annotation/resfam.rules"

# contig classification
include: "rules/contigs_classification/deepvirfinder.rules"
include: "rules/contigs_classification/core_genes.rules"
include: "rules/contigs_classification/plasmids.rules"
include: "rules/contigs_classification/COG_mobilome.rules"
include: "rules/contigs_classification/kaiju.rules"
include: "rules/contigs_classification/cat_bat.rules"
include: "rules/contigs_classification/contigs_statistics.rules"
include: "rules/contigs_classification/combine.rules"

# alignments
include: "rules/genomes_alignment/mummer.rules"