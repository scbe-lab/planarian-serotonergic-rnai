library(data.table)

smed_GOs <- read.delim2(
  file = "/mnt/sda/alberto/projects/smed_rink_gene_annot/predicted.pep.emapper.annotations", skip = 4, header = TRUE, sep = "\t"
)

rosetta <-  read.delim2(
  file = "/mnt/sda/alberto/projects/smed_rink_gene_annot/20231103_Smed_Rink_Simplified_Annotation_Table.tsv", header = TRUE, sep = "\t"
)

source("/mnt/sda/alberto/projects/smed_cisreg/code/r_code/functions/translate_ids.R")

smed_GOs <- smed_GOs[,c(1,10)]
smed_GOs_orig <- smed_GOs
smed_GOs$id <- translate_ids(smed_GOs$X.query,dict = rosetta[c(2,1)])
smed_GOs <- smed_GOs[,c(3,2)]

colnames(smed_GOs) <- c("id","GO")

smed_GOs <- smed_GOs[smed_GOs$GO != "-",]

library(data.table)
setDT(smed_GOs)
smed_GOs <- as.data.frame(smed_GOs[, .(id, GO = unlist(strsplit(as.character(GO), ","))), by = id])[,c(2,3)]


sero_GOs <- read.table("/mnt/sda/elena/sero_manuscript/GO_sero/sero_list_gos.tsv",header = FALSE, sep = "\t")

smed_genes_annotated_with_sero_GOs <- unique(smed_GOs$id[smed_GOs$GO %in% sero_GOs$V1])

length(smed_genes_annotated_with_sero_GOs)

write.table(smed_genes_annotated_with_sero_GOs, "smed_genes_annotated_with_sero_GOs.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)



