# microbiome_shotgun_pipeline


## TODO

### anvi'o test

- [ ] read http://merenlab.org/2016/06/22/anvio-tutorial-v2/
- [ ] see http://merenlab.org/2018/07/09/anvio-snakemake-workflows/
- [ ] gene call: import prodigal gene call
- [ ] kaiju annotation for gene call

### metagenome assembly and binning

- [ ] CAT/BAT taxonomy of contigs
- [ ] see http://merenlab.org/tutorials/infant-gut/#chapter-ii-incorporating-automatic-binning-results
- [ ] binning metabat
- [ ] binning CONCOCT (embedded in anvio)
- [ ] binning groopM
- [ ] binning MyCC
- [ ] binnin BinSanity
- [ ] bin refinment (https://github.com/cmks/DAS_Tool)
- [ ] dereplication (https://github.com/MrOlm/drep)
- [ ] show assembly contiguity on GC vs Coverage plot

### relative abundance

- [X] mOTUs
- [ ] visualization

### AMR, virulence and mobility gene identification

- [ ] VFDB + SwissProt toxins/virulence factors: extract data from fasta and sqlite database
- [X] COG X (Mobilome: prophages, transposons): ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab (119 entries): cdd/COG_mobilome
- [ ] UNIPROT Transposable element [KW-0814] (taxonomy:"Bacteria [2]" keyword:"Transposable element [KW-0814]" AND reviewed:yes) (433 entries)
- [X] CARD (protein)
    - [X] remove variant model, resistance by absence, molecular bypass, ...
- [ ] cluster at 95% identity
- [X] mmseqs2 indexing
- [X] mmseqs2 search
- [X] deepvirfinder
- [X] calculate RPKM (reads per million per sequence length)
- [X] heatmap abundance of VF
- [ ] checkM => get list of core hmm
- [ ] identification of mobile genetic elements
    - [ ] plot contig size of resistance encoding contigs vs non resistance encoding plasmids 
    - [ ] idem with separation by mechanism
    - [ ] BLASTn plasmid database? cutoffs to determine whether is a plasmid or not?
    - [ ] hmm search plasmid database?
    - [ ] extract COG mobilome from COG database, rpsblast COG mobilome only
