library("ggplot2")
library("reshape2")
library("RSQLite")
library("gplots")
library("svglite")
library("dplyr")
library("tibble")
library("gridExtra")

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


# PLOT 1 number of VFs
res <- dbSendQuery(con, "SELECT sample,group_2,group_1,count(*) as n FROM sequence_counts group by sample,group_2,group_1")
table_counts <- dbFetch(res)

p <- ggplot(data=table_counts, aes(x=sample, y=n, fill=group_2))
p <- p + geom_bar(stat="identity")
p <- p+ theme(axis.text.x = element_text(angle = 90))+ facet_grid(. ~ group_1, scales="free")

ggsave(snakemake@output[["plot1"]], p, height=5, width=8)
dev.off()


# PLOT 2
pdf(snakemake@output[["plot2"]])
for (i in unique(sample_table$group_1)) { 
    str <- paste0('select group_2,genus,count(*) as n from (select distinct sample,t1.accession,group_2,genus from sequence_counts t1 inner join uniparc_accession2genus t2 on t1.accession=t2.accession where group_1="', i ,'") A group by group_2,A.genus order by n DESC;')
    res <- dbSendQuery(con, str)
    table_genus <- dbFetch(res)
    print(table_genus)
    table_genus_dcast <- dcast(table_genus, genus~group_2, value.var="n")
    table_genus_dcast[is.na(table_genus_dcast)] <- 0

    dendro <- as.dendrogram(hclust(d = dist(x = as.matrix(table_genus_dcast))))
    order <- order.dendrogram(dendro)

    table_genus$genus <- factor(x = table_genus$genus,
                                levels = table_genus_dcast[,1][order], 
                                ordered = TRUE)

    p <- ggplot(table_genus , aes(group_2, genus)) + geom_tile(aes(fill = n)) + scale_fill_gradient(low = "yellow", high = "steelblue")
    p <- p + theme(axis.text.x = element_text(angle = 90))
    print(p)
    ggsave(paste0("test", i, ".svg"), p, height=9, width=5)
}
dev.off()