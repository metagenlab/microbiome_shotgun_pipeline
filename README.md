# microbiome_shotgun_pipeline

This branch of the microbiome shotgun pipeline includes metagenomic classification tools benchmarking and their database setup.
## Database setup
So far there are 5 tools for which there are rules to generate a database, namely: SURPI, Pathseq, Kraken2, bracken, ganon
and ezVIR.
### Database setup for SURPI
SURPI is ran through a singularity container and requires specific singularity arguments which can be incompatible with other tools using singularity (In our case Pathseq).
To avoid this we run SURPI commands seperately.
To generate the required database files run the command below
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml
 --use-conda --conda-prefix path/to/conda/envs/ 
--use-singularity --singularity-args --app SURPI all_surpi_db
```
### Database setup for the rest of tools
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml 
--use-conda --conda-prefix path/to/conda/envs/ 
--use-singularity --cores 60 all_database_setup 
```

## Getting tools taxonomy

### Getting SURPI outputs
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml
 --use-conda --conda-prefix path/to/conda/envs/ 
--use-singularity --singularity-args --app SURPI 
--bind /data/databases/ --resources mem_mb=1 --cores 100 all_surpi
```
the mem_mb resource controls the number of sample to run in parallel. However, it is advised to run one sample at a time.
### Getting taxonomy for all tools
After getting SURPI outputs, we can now run the rest of the pipeline to generate taxonomy tables for all tools. For this,
run the command below
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml
 --use-conda --conda-prefix path/to/conda/envs/ 
--use-singularity --bind /data/databases/ --resources job_memory_limit=30 --cores 100 all_tax
```
The resource job_memory_limit is the maximum number of RAM in GB allocated for one job. Note that only Pathseq limits its memory usage.
Usually to make sure that jobs run in parallel without exceeding available RAM, we increase the number of cores in the config file.
The number of jobs running in parallel is the number of cores in the snakemake command divided by the number cores specified in the config file.

### Visualizing tool output
We developed scripts to generate heatmaps and other plots to summarize and visualize each tool's output, for both bacteria and viruses.
Define the superkingdom and the taxonomic rank in the config file and then run the command below:
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml
 --use-conda --conda-prefix path/to/conda/envs/ --cores 10 all_tool_out 
```

### Benchmarking according to a gold standard
To compare taxonomy results from each tool to a reference, the user has to create a table with the expected species and their read counts for each sample.
The table has to be created as shown below:
```bash
mkdir pipelinedir/gold_standard
cat pipelinedir/gold_standard/name_of_gold_standard.tsv
sample	taxid	read_counts
AS1_R	42789	2135
AS1_D	443975	87
AS2_D	10372	572
AS2_R	185910	261
AS2_D	687358	6975
AS3_D	10372	1884
```
#### Generating benchmarks for Viruses or Bacteria
Similar to visualizing tool outputs, choose a superkingdom and a taxonomic rank for which the benchmark will be made and run the command below: 
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml
 --use-conda --conda-prefix path/to/conda/envs/ --cores 10 all_benchmark
```


## Coverage plots of reference genomes
The user has to create a directory containing the fasta.gz files of interest:
```bash
mkdir pipelinedir/coverage
cp *.fna.gz > pipelinedir/coverage
```
To generate coverage plots, run the command below:
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml 
--use-conda --conda-prefix path/to/conda/envs/ 
--cores 100 all_coverage -k
```
Here the -k flag ensures that the piepline does not stop when encountering an error. Usually this is due to the fact that 
no read maps to the reference fasta which generates an emtpy BAM file.

## Assembly and read mapping against contigs
If read mapping yields low reference genome coverage, we can assemble reads into contigs, BLASTN contigs against the 
reference genome and use the best mapping contigs to generate coverage plots.

For this, we first assemble and select best matching contigs:
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml 
--use-conda --conda-prefix path/to/conda/envs/ 
--cores 100 all_contig_coverage
```
Then, we can copy the contigs files in contigs-blastn directory to the coverage directory:
```bash
cp reference-accession-name/contigs-blastn/samplename-contigs.fasta ..
gzip samplename-contigs.fasta > samplename-contigs.fna.gz
```
Finally rerun the coverage command to generate the plots:
```bash
snakemake --snakefile path/to/Snakefile --configfile config.yml 
--use-conda --conda-prefix path/to/conda/envs/ 
--cores 100 all_coverage -k
```

