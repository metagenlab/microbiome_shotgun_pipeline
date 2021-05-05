# Title     : Alpha diversity measures
# Objective : correlate alpha diversities between sample types
# Created by: Farid Chaabane
# Created on: 2021-05-03
library('ggplot2')
library('reshape2')
library('plyr')
library('phyloseq')
#load data
readcounts <- read.csv(snakemake@input[[1]], sep=',', header=T)
metadata <- read.csv(snakemake@input[[2]], sep=',', header=T)

#get phyloseq object
otu <- otu_table(readcounts, taxa_are_rows = F)
sam <- sample_data(metadata)
ps <- phyloseq(otu_table=otu, sample_data=sam)

#calculate alpha diversity measures and add metadata
alphadiv <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
alphadiv$samples <- ps@sam_data$sample
alphadiv$types <- ps@sam_data$type
alphadiv$bodysites <- ps@sam_data$bodysite

# reshape table
alphadiv <- melt(alphadiv)
names(alphadiv) <- c('samples','types','bodysites','alpha','value')

# subset sample types, this works only if "type" is a grouping variable with 2 levels (simulated and real for my data)
utypes <- unique(alphadiv$types)
subset1 <- subset(alphadiv, types==utypes[1])
subset2 <- subset(alphadiv, types==utypes[2])

#Get the x and y axis names
type1 <- paste(utypes[1],'_alpha',sep='')
type2 <- paste(utypes[2],'_alpha',sep='')
names(subset1)[names(subset1)=='value'] <- "value1"
names(subset2)[names(subset2)=='value'] <- "value2"

#combine tables
df <- cbind(subset1, subset2)
df <- subset(df, select=which(!duplicated(names(df)))) #remove duplicated columns

# Calculate correlation for each alpha diversity measure
r <- snakemake@params[[1]] #Correlation method
cors <- ddply(df, c("alpha"), summarise, cor=round(cor(value1, value2, method=r), 2),
              x=round(mean(value1)),y=round(max(value1)))  # x and y correspond to correlation plot positions
cors$r2 <- round(cors$cor^2,2) #get r squared
# Plot correlations
png(snakemake@output[[1]], res=400, units = 'in', height = 10, width = 20)
ggplot(df, aes(x=value2, y=value1)) +
  geom_point() +
  xlab(type2) +
  ylab(type1) +
  facet_wrap(~alpha, scales='free', ncol = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text(data=cors, size=3, aes(x=x, y=y, label=paste("r^2=", r2, sep="")))
dev.off()
