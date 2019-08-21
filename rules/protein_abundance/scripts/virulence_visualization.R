library("ggplot2")
library("reshape2")
library("RSQLite")
library("gplots")
library("svglite")
library("dplyr")
library("tibble")

# get table
print(snakemake@input[["database"]])
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

fact <- c("Amsterdam_RAW", "Amsterdam_TREATED", "Evry_RAW", "Evry_TREATED", "Gava_RAW_S1", "Gava_RAW_S2", "Gava_TREATED_S1", "Gava_TREATED_S2", "Prat_RAW", "Prat_TREATED" )

p <- ggplot(data=table_counts, aes(x=sample, y=n, fill=group_2))
p <- p + geom_bar(stat="identity")
p <- p+ theme(axis.text.x = element_text(angle = 90))+ facet_grid(. ~ group_1, scales="free")

ggsave(snakemake@output[["plot1"]], p, height=5, width=8)
