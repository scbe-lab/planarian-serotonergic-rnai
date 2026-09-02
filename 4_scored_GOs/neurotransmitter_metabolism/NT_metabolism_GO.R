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


glut_GOs <- read.table("/mnt/sda/elena/sero_manuscript/GO_sero/glutamate_metabolism.txt",header = FALSE, sep = "\t")
gaba_GOs <- read.table("/mnt/sda/elena/sero_manuscript/GO_sero/GABA_metabolism.txt",header = FALSE, sep = "\t")
choline_GOs <- read.table("/mnt/sda/elena/sero_manuscript/GO_sero/choline_metabolism.txt",header = FALSE, sep = "\t")
amine_GOs <- read.table("/mnt/sda/elena/sero_manuscript/GO_sero/amine_metabolism.txt",header = FALSE, sep = "\t")


smed_genes_annotated_with_glut_GOs <- unique(smed_GOs$id[smed_GOs$GO %in% glut_GOs$V1])
smed_genes_annotated_with_gaba_GOs <- unique(smed_GOs$id[smed_GOs$GO %in% gaba_GOs$V1])
smed_genes_annotated_with_choline_GOs <- unique(smed_GOs$id[smed_GOs$GO %in% choline_GOs$V1])
smed_genes_annotated_with_amine_GOs <- unique(smed_GOs$id[smed_GOs$GO %in% amine_GOs$V1])



length(smed_genes_annotated_with_glut_GOs)
length(smed_genes_annotated_with_gaba_GOs)
length(smed_genes_annotated_with_choline_GOs)
length(smed_genes_annotated_with_amine_GOs)





write.table(smed_genes_annotated_with_glut_GOs, "smed_genes_annotated_with_glut_GOs.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(smed_genes_annotated_with_gaba_GOs, "smed_genes_annotated_with_gaba_GOs.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(smed_genes_annotated_with_choline_GOs, "smed_genes_annotated_with_choline_GOs.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(smed_genes_annotated_with_amine_GOs, "smed_genes_annotated_with_amine_GOs.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)



