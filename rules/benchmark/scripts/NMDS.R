#Redirect console output to log
log <- file(snakemake@log[[1]], open="wt")
       sink(file=log, type=c("output", "message"))

#load libraries
library('vegan')
library('ggplot2')

# load abundance matrix and metadata
readcounts <- read.csv(snakemake@input[[1]], sep=',', header=T)
metadata <- read.csv(snakemake@input[[2]], sep=',', header=T)
#load params
metric <- snakemake@params[[1]]

# Calculate matrix distance
if (metric == 'jaccard') {
  mdist <- vegdist(readcounts, method = metric, binary = T)
} else {
  mdist <- vegdist(readcounts, method = metric, binary = F)
}

# Perform NMDS
metadata$richness <- specnumber(readcounts) #Calculate richness for point size on the plot
nmds <- metaMDS(mdist)
NMDS1 <- nmds$points[,1]
NMDS2 <- nmds$points[,2]
NMDS <- data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, Bodysite = metadata$bodysite,
                   Richness = metadata$richness, Type=metadata$type)
png(snakemake@output[[1]], res=400, units = 'in',height = 5, width = 8)
ggplot(NMDS, aes(x = NMDS1, y = NMDS2, col = Bodysite, label = rownames(NMDS))) +
  geom_point(aes(size = Richness, color = Bodysite, shape = Type)) +
  stat_ellipse(geom = "polygon", aes(group = Bodysite, color = Bodysite, fill = Bodysite), alpha = 0.3) +
  theme_bw() +
  ggtitle(sprintf("NMDS plot, %s distance, %s stress", metric, format(nmds$stress, digits = 4))) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
# PERMANOVA
adonis2(readcounts~type*bodysite, data=metadata , permutations = 999, method=metric)
# Permutest
permutest(betadisper(mdist, group=metadata$type))




