# microbiome_shotgun_pipeline


## TODO

### anvi'o test

- [ ] read http://merenlab.org/2016/06/22/anvio-tutorial-v2/
- [ ] see http://merenlab.org/2018/07/09/anvio-snakemake-workflows/


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

### relative abundance

- [ ] mOTUs
- [ ] visualization

### AMR, virulence and mobility gene identification

- [ ] COG X (Mobilome: prophages, transposons): ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab (119 entries)
- [ ] UNIPROT Transposable element [KW-0814] (taxonomy:"Bacteria [2]" keyword:"Transposable element [KW-0814]" AND reviewed:yes) (433 entries)
- [ ] CARD (protein)
    - [ ] remove variant model, resistance by absence, molecular bypass, ...
- [ ] VFDB + SwissProt toxins/virulence factors
- [ ] cluster at 95% identity
- [ ] mmseqs2 indexing
- [ ] mmseqs2 search
- [ ] calculate RPKM (reads per million per sequence length)
- [ ] heatmap abundance of VF
