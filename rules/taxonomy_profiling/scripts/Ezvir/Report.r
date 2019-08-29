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

# makes 2D plot of genome coverage vs depth using dir of PBC.csv files
# -label names are parsed from the file names
# -dot size scaling based on total genome nucleotides covered
# -colored accoring to virus group
# -also plot genome expect size (ring around dot)
# -opt$color_list = color file
# -opt$j = genome file
# -opt$min_coverage = % genome coverage cutoff (below this value are not plotted)
# -opt$outheader = out file header name (e.g. sample-01)
# -opt$legends = don't print legends
# -opt$label_angle = label angle (0-360)
# -opt$label_size = label size

library('optparse')
library('plyr')
library('TeachingDemos')
suppressMessages(library('zoo'))

### Global ###
RL  = 100         # read length (used to calculate total mapped reads)
RWK = 50          # size of rolling window for depth maximum value
MCL = 400         # minimum covered length (bp), below this, virus not plotted
DBA = 4000        # dot bin size cutoff "A"
DBB = 20000       # dot bin size cutoff "B"
DBC = 40000       # dot bin size cutoff "C"
CTR = 190         # dot color transparency (0 to 255)

PC  = numeric()   # percent covered
DC  = numeric()   # depth of coverage
TNT = numeric()   # total nucleotides covered
TMR = numeric()   # total mapped reads (pileup / read len)
DSC = numeric()   # dot scale factor
DGL = numeric()   # dot genome length ring
GL  = numeric()   # genome lengths of each virus
DCN = character() # dot color name
FN  = character() # file names
LN  = character() # label names (split from file names - short phase-1)
GN  = character() # label names (long genome names for phase-2)
ID  = character() # genome ID

cat('---------------------------------------\n')
cat('ezVIR  Copyright (C) 2014  Tom J. Petty\n')
cat('---------------------------------------\n')
cat('GNU General Public License\n')
cat('<http://www.gnu.org/licenses/>\n\n')

### FUNCTIONS ###
# function for plotting and number padding (to make nice plot axes)
roundUp <- function(x,to=10)
{
  to*(x%/%to + as.logical(x%%to))
}

# adds transparency to a color
# Define transparency with an integer between 0 and 255
# 0 = transparant, 255 = fully visable
# Works with either color and trans a vector of equal length, or one of the two of length 1.
addTrans <- function(color,trans)
{
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

# create options
option_list <- list(
   make_option(c("-c", "--color_list"),
       help="[required] list of colors per virus group"),
   make_option(c("-j", "--genome_list"),
       help="[required] list of virus genomes"),
   make_option(c("-i", "--csvdir"),
       help="[required] input dir of .csv"),
   make_option(c("-o", "--outheader"),
       help="[required] header name for output files"),
   make_option(c("-b", "--blacklist"),
       help="list of virus (shortnames) NOT to plot"),
   make_option(c("-p", "--phase2"),
       help="one virus shortname to plot"),
   make_option(c("-l", "--legends"), action="store_false", default=TRUE,
       help="with this flag, legends are not printed [default=TRUE]"),
   make_option(c("-m", "--min_coverage"), type="integer", default=0,
       help="cutoff for minimum % coverage [default=0]"),
   make_option(c("-s", "--label_size"), type="double", default=1,
       help="size of labels on plot [default=1.0]"),
   make_option(c("-a", "--label_angle"), type="integer", default=0,
       help="angle of data labels on plot"),
   make_option(c("-x", "--xmin"), type="integer", default=0,
       help="min value for x-axis"),
   make_option(c("-z", "--xmax"), type="integer", default=150,
       help="max value for x-axis")
)

# load options from command line
opt <- parse_args(OptionParser(option_list = option_list))

# check that color list exists
if (is.null(opt$color_list)) {
    warning('missing color list [flag: -c]')
    q()
}
# check that genome list exists
if (is.null(opt$genome_list)) {
    warning('missing genome list [flag: -j]')
    q()
}
# check that header name exists
if (is.null(opt$outheader)) {
    warning('missing header name for output files [flag: -o]')
    q()
}
# check that directory of .csv was entered
if (is.null(opt$csvdir)) {
    warning('missing input dir of .csv [flag: -i]')
    q()
}
# check that "-m" is reasonable value [0 to 100]
if (opt$min_coverage > 100) {
    warning(paste('--please select % coverage cutoff in range 0 to 100 [flag: -m]\n',
        '--the value entered was: ',opt$min_coverage,'',sep=""))
    q()
}

### load virus GROUP names and COLORS ###
colors <- read.csv(opt$color_list, header=TRUE, sep=",")
genomes <- read.csv(opt$genome_list, header=TRUE, sep=",")

### load blacklist if supplied ###
blacklist = data.frame()
if (!is.null(opt$blacklist)) {
  blacklist <- read.csv(opt$blacklist, header=FALSE, sep=",")
  # remove any whitespace after each name
  blacklist$V1 <- lapply(blacklist$V1, function(x) gsub("\\s","", x))
}

# give column names to genome list
colnames(genomes) <- c("gid","gname","vgroup","shortname")

### merge GROUP names and COLORS (add corresponding color column) ###
genomes$gcolor <- colors$colorname[match(genomes$vgroup,colors$v_grp)]

### read the directory of .csv, calculate metrics, make plot ###
FN = list.files(normalizePath(opt$csvdir), full.names=T, pattern="*.csv")

# don't use empty files
for (file in FN) {
  FN = FN[file.info(FN)$size > 1]
}

# load in all *csv files
for (i in 1:length(FN)) assign(FN[i], read.csv(FN[i],header=F, sep="\t"))

# gather and make data for plots
for (file in FN) {
  data = get(file)
  colnames(data) <- c("genome","pos","coverage")

  ### get how much of genome is covered ###
  genome_len = max(data$pos) # get length of the genome
  naked = sum(data$coverage == 0) # get genome regions with no reads mapped
  covered = genome_len - naked # bases that are covered
  percent_covered = (covered / genome_len) * 100   # % of genome covered at least once
  mapped_reads = floor(sum(data$coverage) / RL)

  ### get depth as maximum of rolling window of size k
  win_depth <- max(rollmean(data$coverage, RWK))

  # plot only those above the % genome covereage cutoff value AND
  # with more than MCL - minimum covered length nucleotides
  if ((percent_covered >= opt$min_coverage) & covered >= MCL ) {

    ### get short label using this file name ###
    # remove extension (.csv) from genome name
    bn = basename(file)
    gen_name = gsub(".csv","",bn)

    # get row in "genomes" with corresponding short name
    intersect <- match(gen_name, genomes$gname)

    ### catch the case when label names don't match ###
    # using list of names to plot (GENOMES.ezv), can have subsets of all virus tested.
    if (is.na(genomes[intersect,"shortname"])) {
#      print("---Intersect NA DETECTED---")
#      print(gen_name)
    }

    # only update vectors if label is understood
    if (!is.na(genomes[intersect,"shortname"])) {

      ### update vectors ###
      PC  = c(PC, percent_covered)
      DC  = c(DC, win_depth)
      GL  = c(GL, genome_len)
      TNT = c(TNT, covered)
      TMR = c(TMR, mapped_reads)

      # scale the colored dot size accoring to total covered nucleotides
      if (covered < DBA) {
        DSC = c(DSC, 1)
      }
      else if (covered >= DBA & covered < DBB) {
        DSC = c(DSC, 2)
      }
      else if (covered >= DBB & covered < DBC) {
        DSC = c(DSC, 3)
      }
      else {
        DSC = c(DSC, 4)
      }

      # scale the genome length dot size accoring to genome length
      # also scale the label offset (LO) so it's not overlapping bigger dots
      if (genome_len < DBA) {
        DGL = c(DGL, 1)
      }
      else if (genome_len >= DBA & genome_len < DBB) {
        DGL = c(DGL, 2)
      }
      else if (genome_len >= DBB & genome_len < DBC) {
        DGL = c(DGL, 3)
      }
      else {
        DGL = c(DGL, 4)
      }

      # store corresponding color name
      cn = toString(genomes[intersect,"gcolor"])
      # if no color, set to black
      if (is.na(genomes[intersect,"gcolor"])) {
        cn = "black"
      }

      # update color vector for this point
      DCN = c(DCN, cn)

      # store short,long label name and genome ID
      LN = c(LN, toString(genomes[intersect,"shortname"]))
      GN = c(GN, gen_name)
      ID = c(ID, toString(data$genome[1])) # store the genome ID

    } # end if label name != na
  } # end if percent_covered > input value
}

### merge all data ###
all_recs = data.frame(cbind(LN, GN, ID, DCN, DSC, PC, DC, TNT, GL, DGL, TMR))

# if there were no virus with % genome coverage above specified "-m"
# quit nicely, there is nothing to plot!
if (nrow(all_recs) < 1) {
    warning(paste("--No genomes with more than ", opt$min_coverage ," % genome coverage\n",
                  "--Consider lowering the value of '-m'", sep=""))
    q()
}

### turn labels, names into strings ###
all_recs$ID <- as.character(all_recs$ID)
all_recs$LN <- as.character(all_recs$LN)
all_recs$DCN <- as.character(all_recs$DCN)

### convert from factors to real numbers ###
all_recs$PC <- as.numeric(as.character(all_recs$PC))
all_recs$DC <- as.numeric(as.character(all_recs$DC))
all_recs$DSC <- as.numeric(as.character(all_recs$DSC))
all_recs$TNT <- as.numeric(as.character(all_recs$TNT))
all_recs$GL <- as.numeric(as.character(all_recs$GL))
all_recs$DGL <- as.numeric(as.character(all_recs$DGL))
all_recs$TMR <- as.numeric(as.character(all_recs$TMR))

### get unique label names ###
# get a unique list of all short label names
unq_recs <- all_recs[!duplicated(all_recs$LN),]
# store only the unique label names for next step
unq_ln <- unq_recs$LN

### foreach unique (plot best of each family) ###
# store values for plotting (plot values)
pv <- data.frame() # records that will be plotted
bg <- data.frame() # best genomes for all families, even blacklisted


#----------------PHASE-1 ONLY------------------------------------
if(is.null(opt$phase2)) {

    # get the best for each shortname
    # 1) get the best depth of coverage then
    # 2) get the best TNT (if multiple genomes have same depth
    tt = ddply(all_recs,.(LN),function(DF) {res <- DF[which.max(DF$PC),]
                        res[which.max(res$TNT),]})

    # add best genome for each unique label name
    for (uln in unq_ln) {

        # update "all best genomes"
        bg = rbind(bg, subset(tt, LN == uln))

        # don't consider this virus if it's on the blacklist
        if ((nrow(blacklist) > 0) & (uln %in% blacklist$V1)) {
            warning(paste("Did not plot blacklisted virus: ", uln,"", sep=""))
        }

        else {
            pv = rbind(pv, subset(tt, LN == uln))
        } # end "if this virus is not on blacklist"
    }
}
#----------------PHASE-1 ONLY------------------------------------


#----------------PHASE-2 ONLY------------------------------------
# here the short name will be given , take all members of this name
# for each unique label name
if (!is.null(opt$phase2)) {

    # short label name
    sln = opt$phase2

    # make a subset of all records matching given short label name (sln)
    subset <- all_recs[all_recs$LN == sln,]
    print(subset)

    if (nrow(subset) < 1) {
        warning(paste('No records found for Phase-2 name: ', sln,sep=""))
        q()
    }

    else {
        # store this information for plotting
        pv = rbind(pv,subset)
    }
}
#----------------PHASE-2 ONLY------------------------------------


# turn labels, names into strings
colors$v_grp <- as.character(colors$v_grp)
colors$colorname <- as.character(colors$colorname)

# get the color list in alphabetical order (for printing legend on plot)
order.colors <- order(colors$v_grp)
colors <- colors[order.colors, ]

# add "other case" and padding to bottom of list
colors <- rbind(colors, c("Other","black"))
colors <- rbind(colors, c("padding","white"))


### make the plots ###
# = paste("cut-", opt$min_coverage, sep="")
filename = paste(opt$outheader,"p1-plot.pdf", sep = "_")
# if phase-2, write the short name into filename
if (!is.null(opt$phase2)) {
    filename = paste(opt$outheader, covinfo, opt$phase2,"p2-plot.pdf", sep = "_")
}
#if blacklist was used, add 'BL' to filename
if (!is.null(opt$blacklist)) {
    filename = paste(opt$outheader, covinfo, "BL-p1-plot.pdf", sep = "_")
}

pdf(filename)
par(oma=c(1.5,0,0,0))
ylabel = paste("max coverage depth [",RWK, " bp window]", sep = "")

### gray genome size rings ###
# first, plot color dots to show coverage for each point
plot(pv$PC, pv$DC, xlab="", ylab=ylabel,
     abline(h = 0, v = 100, lty = 2, col = "gray60"),
     xlim=c(opt$xmin,opt$xmax) ,ylim=c(0, roundUp(max(pv$DC))),
     xaxt="n", pch=19, cex=pv$DSC, col=addTrans(pv$DCN,CTR))

# second, add the "genome expect size" rings (around dots)
points(pv$PC, pv$DC, pch=1, lwd=1, cex=pv$DGL, col="gray")


### color legends ###
if (opt$legends) {
# legend for dot size
    dotsize_A = paste(" < ", format(DBA, big.mark=",", scientific=FALSE), sep = "")
    dotsize_B = paste(" > ", format(DBA, big.mark=",", scientific=FALSE), sep = "")
    dotsize_C = paste(" > ", format(DBB, big.mark=",", scientific=FALSE), sep = "")
    dotsize_D = paste(" > ", format(DBC, big.mark=",", scientific=FALSE), sep = "")
    legnames = c(dotsize_A, dotsize_B, dotsize_C, dotsize_D)
    legend('topright', legnames, pch=c(19,19,19,19), pt.cex=c(1,2,3,4), col="gray80", text.col="gray80",bty="n", cex=0.7,
         title="Total covered length (bp): ", x.intersp = 1.5, y.intersp = c(1,1,1.15,1.35))

# legend for virus families
    legend('bottomright', colors$v_grp, pch=15, pt.cex=1, col=addTrans(colors$colorname,CTR), text.col=colors$colorname,bty="n", cex=0.7,
         title=" ", title.col="black", x.intersp = 1.5, y.intersp = 1)
}

# adjust x-axis tick marks accordingly
ticks <- c(0,25,50,75,100)
if ((opt$xmin > 0) | (opt$xmax < 150)) {
  ticks <- c(opt$xmin,opt$xmax/4,opt$xmax/2,(opt$xmax/2)+(opt$xmax/4),opt$xmax)
}
axis(side=1, at=ticks)

# data labels (with rotation value as "LR")
if (is.null(opt$phase2)) {
# text(pv$PC, pv$DC, pv$LN, cex=0.6, pos=4, offset=0.9, srt=opt$label_angle, col=addTrans(pv$DCN,CTR), las=2)

  # set up the values for distributing the labels to avoid label overlap on plots
  # this may need to be updated depending on general range of Y-axis
  # - works in most cases!
  ylim=roundUp(max(pv$DC))
  diff=1
  step=1
  if (ylim > 1000) {
    diff=ylim/200
    print(diff)
    step=4
  }
  else if ((ylim > 20) && (ylim < 1000)){
    diff=ylim/50
    step=diff/10
  }
  else {
    diff=0.1
    step=ylim/10
  }
  text(pv$PC, spread.labs(pv$DC, mindiff=diff, maxiter=100000, stepsize=step, min=0, max=ylim), pv$LN, cex=opt$label_size ,pos=4, offset=0.9, srt=opt$label_angle, col=pv$DCN, las=2)

}
# phase-2 use long data point labels
if (!is.null(opt$phase2)) {
#  text(pv$PC, pv$DC, pv$GN, cex=0.6, pos=4, offset=0.9, srt=opt$label_angle, col=addTrans(pv$DCN,CTR), las=2)
  text(pv$PC, pv$DC, pv$GN, cex=1, pos=4, offset=0.9, srt=opt$label_angle, col=pv$DCN, las=2)
}


# main titles and labels
title(opt$outheader, adj=0.3)
subtitle = paste("PHASE-1 | viruses with >", opt$min_coverage, "% genome coverage",  sep = " ")

# change title info for phase-2
if (!is.null(opt$phase2)) {
    subtitle = paste("PHASE-2 |",opt$phase2,"with >", opt$min_coverage, "% genome coverage",  sep = " ")
}

mtext(subtitle,  NORTH <-3, line=0.25, cex=0.7, col="darkgray", font=3, adj=0.2)
mtext(paste(" ",filename, " | ", format(Sys.time(), "%Y-%m-%d %H:%M")),
      cex=0.75, line=0, side=SOUTH<-1, adj=0, outer=TRUE, col="darkgray", font=3)
mtext(side = 1, "% genome coverage", line = 2.5, adj=0.3)

# if necessary, add red vertical cutoff line for opt$min_coverage
if (opt$min_coverage > 0) {
    abline(v = opt$min_coverage, lty = 2, col = "red")
}

# write the plot(s) to PDF
dev.off()

# write values from the points on plot to a table
if (!is.null(opt$phase2)) {
    write.table(pv,file=paste(opt$outheader, covinfo, opt$phase2, "p2-values.txt", sep = "_"), sep = ',', row.names=pv$LN, col.names=NA)
}

if (is.null(opt$phase2)) {
    write.table(pv,file=paste(opt$outheader, covinfo, "p1-values.txt", sep = "_"), sep = ',', row.names=pv$LN, col.names=NA)

    # write BEST genomes names file for use in CC-analysis
    bestinfo = bg[,c("ID","GN","LN", "DC", "DCN")]
    write.table(bestinfo,file="BEST-GENOMES.ezv", sep = ',', row.names=FALSE, quote=FALSE)

    # write all_recs for CC-analysis later
    # don't need to write all the dot scaling info
    # keep the dot color, can be used in plotting
    #  ID  = Genome ID
    #  GN  = Genome Name
    #  LN  = Label Name
    #  GL  = Genome Length
    #  PC  = Percent Coverage
    #  TNT = Total NT (nucleotides) of the genome that are covered (bp)
    #  DC  = Depth of Coverage
    #  DCN = Dot Color Name
    allinfo = all_recs[,c("ID", "GN", "LN", "GL", "PC", "TNT", "DC", "DCN")]
    write.table(allinfo,file="ALL-RESULTS.ezv", sep = ',', row.names=FALSE, quote=FALSE)
}

