library("ggplot2")
library("reshape2")
library("RSQLite")
library("gplots")
library("svglite")
library("dplyr")
library("tibble")
library("gridExtra")

print("PLOTTING 1")
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

print("PLOTTING 2")
# PLOT 1 number of different resistance genes 
res <- dbSendQuery(con, "SELECT sample,group_2,group_1,count(*) as n FROM sequence_counts group by sample,group_2,group_1")
table_counts <- dbFetch(res)

p <- ggplot(data=table_counts, aes(x=sample, y=n, fill=group_2))
p <- p + geom_bar(stat="identity")
p <- p+ theme(axis.text.x = element_text(angle = 90))+ facet_grid(. ~ group_1, scales="free")

ggsave(snakemake@output[[1]], p, height=5, width=8)
dev.off()

print("PLOTTING 3")
# PLOT 2 number of different resistance genes with more than 10 counts
res <- dbSendQuery(con, "SELECT sample,group_2,group_1,count(*) as n from (SELECT sample,group_2,group_1 FROM sequence_counts where n_hits > 9) A group by sample,group_2,group_1")
table_counts <- dbFetch(res)

p <- ggplot(data=table_counts, aes(x=sample, y=n, fill=group_2))
p <- p + geom_bar(stat="identity")
p <- p+ theme(axis.text.x = element_text(angle = 90))+ facet_grid(. ~ group_1, scales="free")

ggsave(snakemake@output[[2]], p, height=5, width=8)
dev.off()

print("PLOTTING 4")
# PLOT 3: distribution of read counts
p <- ggplot(rpkm_table, aes(x = n_hits, fill=group_2))
p <- p + geom_histogram(colour = "white") + facet_grid(group_2 ~ group_1)
ggsave(snakemake@output[[3]], p, height=6, width=8)

print("PLOTTING 4")
# PLOT 4: read counts boxplots
 p <- ggplot(rpkm_table, aes(y = n_hits, x=sample, fill=group_2))
 p <- p + geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) 
 p <- p + theme(axis.text.x = element_text(angle = 90)) + facet_grid(. ~ group_1, scales="free")
 ggsave(snakemake@output[[4]], p, height=6, width=8)

print("PLOTTING 5")
 # PLOT 5 RPKM heatmap all
dendro <- as.dendrogram(hclust(d = dist(x = as.matrix(rpkm_log2_table_dcast))))
order <- order.dendrogram(dendro)

rpkm_table$accession <- factor(x = rpkm_table$accession,
                               levels = rownames(rpkm_log2_table_dcast)[order], 
                               ordered = TRUE)

p <- ggplot(rpkm_table, aes(sample, accession)) + geom_tile(aes(fill = RPKM_log2)) + scale_fill_gradient(low = "white", high = "steelblue")
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + facet_grid( . ~ group_2, scales="free")
ggsave(snakemake@output[[5]], p, height=26, width=8)

print("PLOTTING 6")
# PLOT 6 RPKM heatmap family 
res <- dbSendQuery(con, "select t1.sample,t1.group_1,t1.group_2,t2.AMR_family,SUM(t1.RPKM) as family_sum from sequence_counts t1 inner join accession2aro t2 on t1.accession=t2.protein_accession group by t1.sample,t1.group_1,t1.group_2,t2.AMR_family order by family_sum DESC;
")
AMR_family_RPKM <- dbFetch(res)
AMR_family_RPKM_dcast <- dcast(AMR_family_RPKM, sample~AMR_family, value.var="family_sum")
AMR_family_RPKM_dcast[is.na(AMR_family_RPKM_dcast)] <- 0
# use samw row names and col names
rownames(AMR_family_RPKM_dcast) <- AMR_family_RPKM_dcast[,1]
# transpose the table
AMR_family_RPKM_dcast <- t(AMR_family_RPKM_dcast[,2:length(AMR_family_RPKM_dcast[1,])])
# reorder rows based on rowSum
AMR_family_RPKM_dcast <- AMR_family_RPKM_dcast[order(rowSums(AMR_family_RPKM_dcast),decreasing=F),]
# match index matrix rownames to 
ordered_rows <- rownames(AMR_family_RPKM_dcast)
AMR_family_RPKM$AMR_family <- factor(x = AMR_family_RPKM$AMR_family,
                                    levels = ordered_rows, 
                                    ordered = TRUE)
p <- ggplot(AMR_family_RPKM, aes(sample, AMR_family)) + geom_tile(aes(fill = family_sum)) + scale_fill_gradient(low = "white", high = "steelblue")
p <- p + geom_text(aes(label = round(family_sum, 1)), position = position_dodge(width=0.9),  size=9) 
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + facet_grid( . ~ group_2, scales="free")
ggsave(snakemake@output[[6]], p, height=13, width=20)

print("PLOTTING 7")
# PLOT 7 RPKM heatmap resistance mechanism 
res <- dbSendQuery(con, "select t1.sample,t1.group_1,t1.group_2,t2.resistance_mechanism,SUM(t1.RPKM) as mechanism_sum from sequence_counts t1 inner join accession2aro t2 on t1.accession=t2.protein_accession group by t1.sample,t1.group_1,t1.group_2,t2.resistance_mechanism order by mechanism_sum DESC;;
")
resistance_mechanism_RPKM <- dbFetch(res)
resistance_mechanism_RPKM_dcast <- dcast(resistance_mechanism_RPKM, sample~resistance_mechanism, value.var="mechanism_sum")
resistance_mechanism_RPKM_dcast[is.na(resistance_mechanism_RPKM_dcast)] <- 0
# use samw row names and col names
rownames(resistance_mechanism_RPKM_dcast) <- resistance_mechanism_RPKM_dcast[,1]
# transpose the table
resistance_mechanism_RPKM_dcast <- t(resistance_mechanism_RPKM_dcast[,2:length(resistance_mechanism_RPKM_dcast[1,])])
# reorder rows based on rowSum
resistance_mechanism_RPKM_dcast <- resistance_mechanism_RPKM_dcast[order(rowSums(resistance_mechanism_RPKM_dcast),decreasing=F),]
# match index matrix rownames to 
ordered_rows <- rownames(resistance_mechanism_RPKM_dcast)
resistance_mechanism_RPKM$resistance_mechanism <- factor(x = resistance_mechanism_RPKM$resistance_mechanism,
                                    levels = ordered_rows, 
                                    ordered = TRUE)
p <- ggplot(resistance_mechanism_RPKM, aes(sample, resistance_mechanism)) + geom_tile(aes(fill = mechanism_sum)) + scale_fill_gradient(low = "white", high = "steelblue")
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + facet_grid( . ~ group_2, scales="free")
ggsave(snakemake@output[[7]], p, height=12, width=20)

print("PLOTTING 8")
# PLOT 8 RPKM heatmap resistance mechanism 
res <- dbSendQuery(con, "select t1.sample,t1.group_1,t1.group_2,drug_class, sum(t1.RPKM) as drug_class_RPKM from sequence_counts t1 inner join accession2aro t2 on t1.accession=t2.protein_accession inner join aro_accession2drug_class t3 on t2.aro_accession=t3.aro_accession group by t1.sample,t1.group_1,t1.group_2,drug_class;")
drug_class_RPKM <- dbFetch(res)
drug_class_RPKM_dcast <- dcast(drug_class_RPKM, sample~drug_class, value.var="drug_class_RPKM")
drug_class_RPKM_dcast[is.na(drug_class_RPKM_dcast)] <- 0
# use samw row names and col names
rownames(drug_class_RPKM_dcast) <- drug_class_RPKM_dcast[,1]
# transpose the table
drug_class_RPKM_dcast <- t(drug_class_RPKM_dcast[,2:length(drug_class_RPKM_dcast[1,])])
# reorder rows based on rowSum
drug_class_RPKM_dcast <- drug_class_RPKM_dcast[order(rowSums(drug_class_RPKM_dcast),decreasing=F),]
# match index matrix rownames to 
ordered_rows <- rownames(drug_class_RPKM_dcast)
drug_class_RPKM$drug_class <- factor(x = drug_class_RPKM$drug_class,
                                    levels = ordered_rows, 
                                    ordered = TRUE)
p <- ggplot(drug_class_RPKM, aes(sample, drug_class)) + geom_tile(aes(fill = drug_class_RPKM)) + scale_fill_gradient(low = "white", high = "steelblue")
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + facet_grid( . ~ group_2, scales="free")
ggsave(snakemake@output[[8]], p, height=12, width=20)