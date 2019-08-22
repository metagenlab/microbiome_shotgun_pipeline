library("ggplot2")
library("reshape2")
library("RSQLite")
library("gplots")
library("svglite")
library("dplyr")
library("tibble")
library("gridExtra")

print("PLOTTING-------------")
# get sample table
sample_table <- read.csv(snakemake@params[["sample_table"]], sep='\t')

# get table
con <- dbConnect(RSQLite::SQLite(), snakemake@input[["database"]])
res <- dbSendQuery(con, "SELECT * FROM sequence_counts")
rpkm_table <- dbFetch(res)

# claculate RPKM log2
rpkm_table$RPKM_log2 <- log2(rpkm_table$RPKM + 1)

# prepare matrix of RPKM_log2
rpkm_log2_table_dcast <- dcast(rpkm_table, sample~accession, value.var="RPKM_log2")
# replace na by 0
rpkm_log2_table_dcast[is.na(rpkm_log2_table_dcast)] <- log2(1)
# use samw row names and col names
rownames(rpkm_log2_table_dcast) <- rpkm_log2_table_dcast[,1]
# transpose the table
rpkm_log2_table_dcast <- t(rpkm_log2_table_dcast[,2:length(rpkm_log2_table_dcast[1,])])


# PLOT 1 number of different resistance genes 
res <- dbSendQuery(con, "SELECT sample,group_2,group_1,count(*) as n FROM sequence_counts group by sample,group_2,group_1")
table_counts <- dbFetch(res)

p <- ggplot(data=table_counts, aes(x=sample, y=n, fill=group_2))
p <- p + geom_bar(stat="identity")
p <- p+ theme(axis.text.x = element_text(angle = 90))+ facet_grid(. ~ group_1, scales="free")

ggsave(snakemake@output[[1]], p, height=5, width=8)
dev.off()

# PLOT 2 number of different resistance genes with more than 10 counts
res <- dbSendQuery(con, "SELECT sample,group_2,group_1,count(*) as n from (SELECT sample,group_2,group_1 FROM sequence_counts where n_hits > 9) A group by sample,group_2,group_1")
table_counts <- dbFetch(res)

p <- ggplot(data=table_counts, aes(x=sample, y=n, fill=group_2))
p <- p + geom_bar(stat="identity")
p <- p+ theme(axis.text.x = element_text(angle = 90))+ facet_grid(. ~ group_1, scales="free")

ggsave(snakemake@output[[2]], p, height=5, width=8)
dev.off()

# PLOT 3: distribution of read counts
p <- ggplot(rpkm_table, aes(x = n_hits, fill=group_2))
p <- p + geom_histogram(colour = "white") + facet_grid(group_2 ~ group_1)
ggsave(snakemake@output[[3]], p, height=6, width=8)

# PLOT 4: read counts boxplots
 p <- ggplot(rpkm_table, aes(y = n_hits, x=sample, fill=group_2))
 p <- p + geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) 
 p <- p + theme(axis.text.x = element_text(angle = 90)) + facet_grid(. ~ group_1, scales="free")
 ggsave(snakemake@output[[4]], p, height=6, width=8)


 # PLOT 5 RPKM heatmap all
dendro <- as.dendrogram(hclust(d = dist(x = as.matrix(rpkm_log2_table_dcast))))
order <- order.dendrogram(dendro)

rpkm_table$accession <- factor(x = rpkm_table$accession,
                               levels = rownames(rpkm_log2_table_dcast)[order], 
                               ordered = TRUE)

p <- ggplot(rpkm_table, aes(sample, accession)) + geom_tile(aes(fill = RPKM_log2)) + scale_fill_gradient(low = "white", high = "steelblue")
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + facet_grid( . ~ group_2, scales="free")
ggsave(snakemake@output[[5]], p, height=6, width=8)


# PLOT 6 RPKM heatmap family 
res <- dbSendQuery(con, "select t1.sample,t1.group_1,t1.group_2,t2.AMR_family,SUM(t1.RPKM) as family_sum from sequence_counts t1 inner join accession2aro t2 on t1.accession=t2.protein_accession group by t1.sample,t1.group_1,t1.group_2,t2.AMR_family order by family_sum DESC;
")
AMR_family_RPKM <- dbFetch(res)
AMR_family_RPKM_dcast <- dcast(AMR_family_RPKM, sample~AMR_family, value.var="family_sum")
AMR_family_RPKM_dcast[is.na(AMR_family_RPKM_dcast)] <- 0
# use samw row names and col names
rownames(AMR_family_RPKM_dcast) <- AMR_family_RPKM_dcast[,1]
print(AMR_family_RPKM_dcast)
# transpose the table
AMR_family_RPKM_dcast <- t(AMR_family_RPKM_dcast[,2:length(AMR_family_RPKM_dcast[1,])])
print(AMR_family_RPKM_dcast)
# reorder rows based on rowSum
AMR_family_RPKM_dcast <- AMR_family_RPKM_dcast[order(rowSums(AMR_family_RPKM_dcast),decreasing=T),]
# match index matrix rownames to 
ordered_rows <- rownames(AMR_family_RPKM_dcast)
AMR_family_RPKM$AMR_family <- factor(x = AMR_family_RPKM$AMR_family,
                                    levels = ordered_rows, 
                                    ordered = TRUE)
print("ok2")
p <- ggplot(AMR_family_RPKM, aes(sample, AMR_family)) + geom_tile(aes(fill = family_sum)) + scale_fill_gradient(low = "white", high = "steelblue")
print("ok3")
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + facet_grid( . ~ group_2, scales="free")
ggsave(snakemake@output[[6]], p, height=12, width=20)
