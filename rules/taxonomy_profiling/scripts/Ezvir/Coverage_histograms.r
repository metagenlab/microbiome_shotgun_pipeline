#!/usr/bin/Rscript

# =======================================
# ezVIR  Copyright (C) 2014  Tom J. Petty
# =======================================
#
# This file is part of ezVIR-v0.1.
#
#    ezVIR is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ezVIR is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ezVIR.  If not, see <http://www.gnu.org/licenses/>.
#    Please refer to the file "COPYING"

# makes 2D histogram of depth of coverage along genome length using PBC.csv files
# label names are parsed from the file name

library(ggplot2)
library(plyr)
library(optparse)

cat('---------------------------------------\n')
cat('ezVIR  Copyright (C) 2014  Tom J. Petty\n')
cat('---------------------------------------\n')
cat('GNU General Public License\n')
cat('<http://www.gnu.org/licenses/>\n\n')

args <- commandArgs(trailingOnly = TRUE)
PC = numeric() # percent covered
DC = numeric() # depth of coverage
FN = character() # file names
ALPHA = 0.6 # transparenecy for fill color

# create option flags
option_list <- list(
   make_option(c("-d", "--dir"),
       help="[required] directory of .csv mapping results"),
   make_option(c("-b", "--best_genomes"),
       help="[required ] BEST-GENOMES.ezv file"),
   make_option(c("-s", "--shortname"),
       help="[required ] one virus shortname to plot"),
   make_option(c("-n", "--table_names"),
       help="[required ] table regrouping virus names"),
    make_option(c("-g", "--genome_length"),
       help="[required] table regrouping genome lengths")
)

opt <- parse_args(OptionParser(option_list = option_list))


# check that directory is okay
if (is.null(opt$dir)) {
    warning('missing dir of mapping results (.csv files) [flag: -d]')
    q()
}
# check that BEST-GENOMES.ezv file was given
if (is.null(opt$best_genomes)) {
    warning('missing BEST-GENOMES.ezv [flag: -b]')
    q()
}

# check that short name was given
if (is.null(opt$shortname)) {
    warning('missing genome shortname name [flag: -s]')
    q()
}

# read the directory of CCA results
fullpath = normalizePath(opt$dir)
FN = list.files(path=fullpath ,pattern="*.csv", full.names=TRUE)

# don't use empty files
for (file in FN) {
  FN = FN[file.info(FN)$size > 1]
}

# load in all mapping .csv files
for (i in 1:length(FN)) assign(FN[i], read.csv(FN[i],header=F, sep="\t"))

# load the "BEST-GENOMES.ezv" file, output from Phase-1
bg_fp = normalizePath(opt$best_genomes)
best_genomes =  read.csv(bg_fp, header=T, sep=",")
common <- subset(best_genomes, LN == opt$shortname)

BEST = toString(common$ID) # sample gid of strongest signal in family
bname = toString(common$GN) #sample name of strongest signal in family
bval = as.integer(common$DC)
bcol = toString(common$DCN) # color name for this virus family


names_tab=read.csv(opt$table_names,sep='\t',header=T)
genome_lengths=read.csv(opt$genome_length,header=TRUE, sep="\t")
# make outfile
pdf(paste(opt$shortname, "coverage-histogram.pdf", sep = "_"), onefile=TRUE, width=10, height=5)

# create histograms for each mapping
for (file in FN) {
  data = read.table(file, header=F, sep="\t")
  colnames(data) <- c("genome","pos","coverage")
  bn=basename(file)
  genname=gsub(".csv","",bn)
  genome_len = genome_lengths[genome_lengths$X==genname,"genome_length"]# get length of the genome
  depth <- sum(data[,"coverage"])/genome_len
  this_genome = toString(data$genome[1])

  ### get how much of genome is covered ###
   # get genome regions with no reads mapped
  covered = length(data$coverage) # get bases that are covered
  percent_covered = (covered / genome_len) * 100   # % of genome covered at least once

  ### get date and time and format nicely ###
  t1 = format(Sys.Date(), format="%B %d, %Y")
  t2 = format(Sys.time(), format="%I:%M %p")

  ### make plot title from this file name ###
  # get rid of linking characters for nicer titles
  name=as.character(names_tab[names_tab$X==genname,"genome_names"])
  cleanname = gsub("_", " ", name)


  plot.title = paste("reads mapped to:", cleanname, sep=" ")
  plot.subtitle = paste("specimen:", opt$ref, ", genome: ",this_genome,",", t1, "at", t2, sep=" ")
  tit = bquote(atop(.(plot.title), atop(.(plot.subtitle), "")))

  # axis labels
  x.upper = paste("coverage: ", as.integer(percent_covered), "%", sep="")
  x.lower = paste("genome length: ", genome_len,
                  " bp | average depth: ", as.integer(depth),
                  " | tot covered bp: ", as.integer(covered), sep="")
  xinfo = bquote(atop(.(x.upper), atop(italic(.(x.lower)), "")))

  data$coverage <- replace(data$coverage, length(data$coverage), 0)

  print(ggplot(data) + geom_area(aes(data$pos,data$coverage), color="black", fill=bcol, alpha=ALPHA)
        + labs(title=tit)
        + theme(plot.title=element_text(size=rel(1.2), colour="black"))
        + xlab(xinfo) + theme(axis.title.x=element_text(size=rel(1.2)))
        + ylab("depth") + theme(axis.title.y=element_text(size=rel(1.2)))
        )
}

dev.off()
