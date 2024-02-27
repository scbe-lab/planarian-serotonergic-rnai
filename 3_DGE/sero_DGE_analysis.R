# SERO DGE ANALYSIS

# Setup
dir <- "/mnt/sda/elena/R_analysis/sero_dge_202310/sero_DGE_per_cluser/"
setwd(dir)

source("/mnt/sda/alberto/colabos/sizes/code/R_code/r_functions/sourcefolder.R")
sourceFolder("/mnt/sda/alberto/colabos/sizes/code/R_code/r_functions/",recursive = TRUE)

library(Matrix)
library(topGO)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggvenn)

# Load necessary data
sero_Idents <- read.delim2("06092023_smed_sero_identities.csv", sep = "," , header = TRUE)
sero_genes <- read.delim2("06092023_smed_sero_genes.csv", sep = "," , header = TRUE)
sero_leiden_col <- read.delim2("leiden_3_colors.csv", sep = ",")
sero_Idents$color <- translate_ids(x=sero_Idents$leiden_3,dict = sero_leiden_col)

sero_ctypes <- unique(sero_Idents[,c("leiden_3","leiden_3_names","broad_names","color")])
sero_ctypes <- merge(
  sero_ctypes,
  as.data.frame(table(sero_Idents$leiden_3)),
  by.x = 1,
  by.y = 1,
  all.x = TRUE
)

rosetta <- read.delim2("/mnt/sda/alberto/projects/smed_rink_gene_annot/20230915_Smed_Rink_Simplified_Annotation_Table.tsv",header = TRUE,sep="\t")

# Load the matrix
sero_X <- readMM(file = "matrix.mtx")
colnames(sero_X) <-  read.table("barcodes.tsv")[,1]
rownames(sero_X) <-  read.table("features.tsv")[,1]
sero_X <- sero_X[rownames(sero_X) %in% sero_genes$X,colnames(sero_X) %in% sero_Idents$X]

# Create table of counts per condition
rep_group <- data.frame(
  rep = c("B1_T1", "B1_T2", "B2_T1", "B2_T2"),
  group = c("A","B","C","D")
)

sero_Idents$group <- 
  translate_ids(
    x = paste(sero_Idents$biological_rep,sero_Idents$technical_rep, sep = "_"),
    dict = rep_group
  )

sero_pseudobulk_cond_reps <- pseudobulk_cond_rep(
  x = sero_X,
  identities = sero_Idents$leiden_3_names,
  conditions = sero_Idents$sample,
  replicates = sero_Idents$group
)


# Clean sampletable

sero_sampletable <- sero_pseudobulk_cond_reps$sampletable

sero_sampletable_lhx = sero_sampletable[sero_sampletable$condition %in% c("gfp(RNAi)", "lhx1/5-1(RNAi)"), ] # Change based on your conditions
sero_sampletable_lhx <- clean_sampletable(sero_sampletable_lhx)

sero_sampletable_pitx = sero_sampletable[sero_sampletable$condition %in% c("gfp(RNAi)", "pitx(RNAi)"), ] # Change based on your conditions
sero_sampletable_pitx <- clean_sampletable(sero_sampletable_pitx)

sero_matrix_lhx <- 
  sero_pseudobulk_cond_reps$matrix[
    rownames(sero_pseudobulk_cond_reps$matrix) %in% rosetta$gene[rosetta$gene_type == "hconf"],
    colnames(sero_pseudobulk_cond_reps$matrix) %in% sero_sampletable_lhx$sample
  ]

sero_matrix_pitx <- 
  sero_pseudobulk_cond_reps$matrix[
    rownames(sero_pseudobulk_cond_reps$matrix) %in% rosetta$gene[rosetta$gene_type == "hconf"],
    colnames(sero_pseudobulk_cond_reps$matrix) %in% sero_sampletable_pitx$sample
  ]


# SINGLE CELL ANALYSIS OVER ALL CELL TYPES
sero_DGE_all_lhx <- list()
for(i in unique(sero_sampletable_lhx$ctype)){
  print(paste0("starting with cell type ",i))
  sero_DGE_all_lhx[[i]] <-
    deseq_pseudobulk(
      count_matrix = sero_matrix_lhx,
      samples_info = sero_sampletable_lhx[,-2],
      celltype = i,
      filter_by = "pvalue", p_threshold = 0.05,
      contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
      plot_results = FALSE, min_passing_samples = 2,
      min_counts_per_sample = 1,
      keep_dubious = FALSE
    )
  print(paste0("done cell type ",i))
}

sero_DGE_all_pitx <- list()
for(i in unique(sero_sampletable_pitx$ctype)){
  print(paste0("starting with cell type ",i))
  sero_DGE_all_pitx[[i]] <-
    deseq_pseudobulk(
      count_matrix = sero_matrix_pitx,
      samples_info = sero_sampletable_pitx[,-2],
      celltype = i,
      filter_by = "pvalue", p_threshold = 0.05,
      contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
      plot_results = FALSE, min_passing_samples = 2,
      min_counts_per_sample = 1,
      keep_dubious = FALSE
    )
  print(paste0("done cell type ",i))
}

# METAPLOT DATA FRAME
sero_diffreg_lhx <- data.frame(
  num_diff = sapply(
    sero_DGE_all_lhx,
    function(x){
      a = x$diffgenes
      if(is.na(a[1])){
        b = 0
      } else{
        b = length(a)
      }
      return(b)
    }
  )
)

sero_diffreg_lhx <- 
  merge(
    sero_diffreg_lhx,
    sero_ctypes,
    by.x = 0,
    by.y = "leiden_3_names",
    all.x = TRUE
  )
colnames(sero_diffreg_lhx)[1] <- "ctype"
sero_diffreg_lhx <- sero_diffreg_lhx[-c(grep("unannotated",sero_diffreg_lhx$ctype)),]

sero_diffreg_pitx <- data.frame(
  num_diff = sapply(
    sero_DGE_all_pitx,
    function(x){
      a = x$diffgenes
      if(is.na(a[1])){
        b = 0
      } else{
        b = length(a)
      }
      return(b)
    }
  )
)

sero_diffreg_pitx <- 
  merge(
    sero_diffreg_pitx,
    sero_ctypes,
    by.x = 0,
    by.y = "leiden_3_names",
    all.x = TRUE
  )
colnames(sero_diffreg_pitx)[1] <- "ctype"
sero_diffreg_pitx <- sero_diffreg_pitx[-c(grep("unannotated",sero_diffreg_pitx$ctype)),]


##REGRESSION
sero_lhx_lm <- lm(sero_diffreg_lhx$num_diff~sero_diffreg_lhx$Freq)
sero_lhx_lm_top <- sero_diffreg_lhx$ctype[sero_lhx_lm$residuals > quantile(sero_lhx_lm$residuals,0.9)]
#sero_diffreg_lhx$is_top <- ifelse(sero_diffreg_lhx$ctype %in% sero_lhx_lm_top, 1,0)

sero_pitx_lm <- lm(sero_diffreg_pitx$num_diff~sero_diffreg_pitx$Freq)
sero_pitx_lm_top <- sero_diffreg_pitx$ctype[sero_pitx_lm$residuals > quantile(sero_pitx_lm$residuals,0.9)]
#sero_diffreg_pitx$is_top <- ifelse(sero_diffreg_pitx$ctype %in% sero_pitx_lm_top, 1,0)


##META PLOT NUM DIFF GENES AND CLUSTER SIZE
lhx_metaplot <- ggplot(sero_diffreg_lhx, aes(x = log(Freq), y = num_diff, label = ctype, color = color)) +
  geom_point() +
  labs(
    title = "Sero cluster size vs # DiffReg lhx1/5-1(RNAi)",
    x = "cluster size (log n cells)",
    y = "no. diffreg genes"
  ) +
  scale_color_identity()+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  geom_text_repel(cex=3) + labs(title = "Cluster size vs # DiffReg lhx1/5-1(RNAi)")

lhx_metaplot
sero_lhx_lm_top

par(mfrow = c(2,2))
plot(sero_diffreg_lhx$num_diff~sero_diffreg_lhx$Freq, bg = sero_diffreg_lhx$color, pch = 21)
plot(sero_lhx_lm$residuals, bg = sero_diffreg_lhx$color, pch = 21)
hist(sero_lhx_lm$residuals, breaks = 20)
abline(v=quantile(sero_lhx_lm$residuals,0.9), lwd = 2, col = "red")
plot(density(sero_lhx_lm$residuals))
abline(v=quantile(sero_lhx_lm$residuals,0.9), lwd = 2, col = "red")
par(mfrow = c(1,1))


pitx_metaplot <- ggplot(sero_diffreg_pitx, aes(x = log(Freq), y = num_diff, label = ctype, color = color)) +
  geom_point() +
  labs(
    title = "Sero cluster size vs # DiffReg pitx(RNAi)",
    x = "cluster size (log n cells)",
    y = "no. diffreg genes"
  ) +
  scale_color_identity()+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  geom_text_repel(cex=3) + labs(title = "Cluster size vs # DiffReg pitx(RNAi)")

pitx_metaplot
sero_pitx_lm_top

par(mfrow = c(2,2))
plot(sero_diffreg_pitx$num_diff~sero_diffreg_pitx$Freq, bg = sero_diffreg_pitx$color, pch = 21)
plot(sero_pitx_lm$residuals, bg = sero_diffreg_pitx$color, pch = 21)
hist(sero_pitx_lm$residuals, breaks = 20)
abline(v=quantile(sero_pitx_lm$residuals,0.9), lwd = 2, col = "red")
plot(density(sero_pitx_lm$residuals))
abline(v=quantile(sero_pitx_lm$residuals,0.9), lwd = 2, col = "red")
par(mfrow = c(1,1))


# Modify p1 and p2 to show expression of a specific gene on the x-axis and log-transformed cluster size as circle size

pseudobulk_matrix <- pseudobulk(x = as.matrix(sero_X), identities = sero_Idents$leiden_3_names)

lhx <- "h1SMcG0003793"
pitx <- "h1SMcG0012776"

lhx_exp <- data.frame(
  ctype = names(pseudobulk_matrix[lhx,]),
  lhx_counts = pseudobulk_matrix[lhx,]
)

sero_diffreg_lhx <- 
  merge(
    sero_diffreg_lhx,
    lhx_exp,
    by.x = 1,
    by.y = "ctype",
    all.x = TRUE
  )


pitx_exp <- data.frame(
  ctype = names(pseudobulk_matrix[pitx,]),
  pitx_counts = pseudobulk_matrix[pitx,]
)

sero_diffreg_pitx <- 
  merge(
    sero_diffreg_pitx,
    pitx_exp,
    by.x = 1,
    by.y = "ctype",
    all.x = TRUE
  )

lhx_metaplot_lhx_expr <- ggplot(sero_diffreg_lhx, aes(x = lhx_counts, y = num_diff, label = ctype, size = Freq, color = color)) +
  geom_point() +
  labs(
    title = "Sero cluster size vs Expression of lhx1/5-1",
    x = "Expression of lhx1/5-1 (n counts)",
    y = "No. DiffReg Genes"
  ) +
  scale_color_identity() +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75)) +
  geom_text_repel(cex = 3) +
  labs(title = "Cluster size vs Expression of lhx1/5-1")

pitx_metaplot_pitx_expr <- ggplot(sero_diffreg_pitx, aes(x = pitx_counts, y = num_diff, label = ctype, size = Freq, color = color)) +
  geom_point() +
  labs(
    title = "Sero cluster size vs Expression of pitx",
    x = "Expression of pitx (n counts)",
    y = "No. DiffReg Genes"
  ) +
  scale_color_identity() +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75)) +
  geom_text_repel(cex = 3) +
  labs(title = "Cluster size vs Expression of pitx")

lhx_metaplot_lhx_expr
pitx_metaplot_pitx_expr


##VOLCANO PLOTS and HEATMAPS
sero_lhx_lm_top
sero_pitx_lm_top

sero_DGE_idvm_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "intestinal and DV muscle",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
idvm_lhx <- sero_DGE_idvm_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "intestinal and DV muscle DEGs lhx1/5-1(RNAi)")
idvm_lhx


sero_DGE_idvm_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "intestinal and DV muscle",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
idvm_pitx <- sero_DGE_idvm_pitx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "intestinal and DV muscle DEGs pitx(RNAi)")
idvm_pitx


sero_DGE_lep_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "late epidermal progenitors",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
lep_lhx <- sero_DGE_lep_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "late epidermal progenitors DEGs lhx1/5-1(RNAi)")
lep_lhx


sero_DGE_lep_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "late epidermal progenitors",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
lep_pitx <- sero_DGE_lep_pitx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "late epidermal progenitors DEGs pitx(RNAi)")
lep_pitx


sero_DGE_pgrn1_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "pgrn+ parenchymal cells 1",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
pgrn1_lhx <- sero_DGE_pgrn1_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "pgrn+ parenchymal cells 1 DEGs lhx1/5-1(RNAi)")
pgrn1_lhx


sero_DGE_pgrn1_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "pgrn+ parenchymal cells 1",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
pgrn1_pitx <- sero_DGE_pgrn1_pitx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "pgrn+ parenchymal cells 1 DEGs pitx(RNAi)")
pgrn1_pitx


sero_DGE_tph_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "tph+ neurons",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
tph_pitx <- sero_DGE_tph_pitx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "gray", fill=NA, linewidth = 0.75))+
  ggtitle(label = "tph+ neurons DEGs pitx(RNAi)")
tph_pitx


sero_DGE_tph_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "tph+ neurons",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
tph_lhx <- sero_DGE_tph_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "tph+ neurons DEGs lhx1/5-1(RNAi)")
tph_lhx


sero_DGE_otf1_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "otf+ neurons 1",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
otf1_lhx <- sero_DGE_otf1_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "otf+ neurons 1 DEGs lhx1/5-1(RNAi)")
otf1_lhx


sero_DGE_otf2_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "otf+ neurons 2",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
otf2_lhx <- sero_DGE_otf2_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "otf+ neurons 2 DEGs lhx1/5-1(RNAi)")
otf2_lhx


sero_DGE_idvm_lhx$heatmap
sero_DGE_idvm_pitx$heatmap

sero_DGE_lep_lhx$heatmap
sero_DGE_lep_pitx$heatmap


sero_DGE_pgrn1_lhx$heatmap
sero_DGE_pgrn1_pitx$heatmap

sero_DGE_tph_pitx$heatmap
sero_DGE_tph_lhx$heatmap

sero_DGE_otf1_lhx$heatmap
sero_DGE_otf2_lhx$heatmap


# GO enrichment analysis

source("/mnt/sda/alberto/colabos/sizes/code/R_code/r_functions/sourcefolder.R")
sourceFolder("/mnt/sda/alberto/colabos/sizes/code/R_code/r_functions/",recursive = TRUE)
source("/mnt/sda/alberto/projects/smed_cisreg/code/r_code/functions/topGO_wrapper.R")

smed_id_GO <- readMappings("/mnt/sda/alberto/projects/smed_cisreg/outputs/gene_annotation/smed_GOs.tsv")

sero_DGE_idvm_lhx$res
sero_DGE_idvm_pitx$res
sero_DGE_lep_lhx$res
sero_DGE_lep_pitx$res
sero_DGE_pgrn1_lhx$res
sero_DGE_pgrn1_pitx$res
sero_DGE_tph_pitx$res
sero_DGE_tph_lhx$res
sero_DGE_otf1_lhx$res
sero_DGE_otf2_lhx$res

significant_idvm_lhx <- subset(sero_DGE_idvm_lhx$res, pvalue < 0.05)
significant_idvm_pitx <- subset(sero_DGE_idvm_pitx$res, pvalue < 0.05)
significant_lep_lhx <- subset(sero_DGE_lep_lhx$res, pvalue < 0.05)
significant_lep_pitx <- subset(sero_DGE_lep_pitx$res, pvalue < 0.05)
significant_pgrn1_lhx <- subset(sero_DGE_pgrn1_lhx$res, pvalue < 0.05)
significant_pgrn1_pitx <- subset(sero_DGE_pgrn1_pitx$res, pvalue < 0.05)
significant_tph_lhx <- subset(sero_DGE_tph_lhx$res, pvalue < 0.05)
significant_tph_pitx <- subset(sero_DGE_tph_pitx$res, pvalue < 0.05)
significant_otf1_lhx <- subset(sero_DGE_otf1_lhx$res, pvalue < 0.05)
significant_otf2_lhx <- subset(sero_DGE_otf2_lhx$res, pvalue < 0.05)

significant_genes_stat_idvm_lhx <- data.frame(significant_idvm_lhx)
significant_genes_stat_idvm_pitx <- data.frame(significant_idvm_pitx)
significant_genes_stat_lep_lhx <- data.frame(significant_lep_lhx)
significant_genes_stat_lep_pitx <- data.frame(significant_lep_pitx)
significant_genes_stat_pgrn1_lhx <- data.frame(significant_pgrn1_lhx)
significant_genes_stat_pgrn1_pitx <- data.frame(significant_pgrn1_pitx)
significant_genes_stat_tph_lhx <- data.frame(significant_tph_lhx)
significant_genes_stat_tph_pitx <- data.frame(significant_tph_pitx)
significant_genes_stat_otf1_lhx <- data.frame(significant_otf1_lhx)
significant_genes_stat_otf2_lhx <- data.frame(significant_otf2_lhx)


list_diffregs_forGO_idvm_lhx <- list(
  down = rownames(significant_genes_stat_idvm_lhx[significant_genes_stat_idvm_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_idvm_lhx[significant_genes_stat_idvm_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_idvm_pitx <- list(
  down = rownames(significant_genes_stat_idvm_pitx[significant_genes_stat_idvm_pitx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_idvm_pitx[significant_genes_stat_idvm_pitx$log2FoldChange > 0,])
)
list_diffregs_forGO_lep_lhx <- list(
  down = rownames(significant_genes_stat_lep_lhx[significant_genes_stat_lep_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_lep_lhx[significant_genes_stat_lep_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_lep_pitx <- list(
  down = rownames(significant_genes_stat_lep_pitx[significant_genes_stat_lep_pitx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_lep_pitx[significant_genes_stat_lep_pitx$log2FoldChange > 0,])
)
list_diffregs_forGO_pgrn1_lhx <- list(
  down = rownames(significant_genes_stat_pgrn1_lhx[significant_genes_stat_pgrn1_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_pgrn1_lhx[significant_genes_stat_pgrn1_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_pgrn1_pitx <- list(
  down = rownames(significant_genes_stat_pgrn1_pitx[significant_genes_stat_pgrn1_pitx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_pgrn1_pitx[significant_genes_stat_pgrn1_pitx$log2FoldChange > 0,])
)
list_diffregs_forGO_tph_lhx <- list(
  down = rownames(significant_genes_stat_tph_lhx[significant_genes_stat_tph_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_tph_lhx[significant_genes_stat_tph_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_tph_pitx <- list(
  down = rownames(significant_genes_stat_tph_pitx[significant_genes_stat_tph_pitx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_tph_pitx[significant_genes_stat_tph_pitx$log2FoldChange > 0,])
)
list_diffregs_forGO_otf1_lhx <- list(
  down = rownames(significant_genes_stat_otf1_lhx[significant_genes_stat_otf1_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_otf1_lhx[significant_genes_stat_otf1_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_otf2_lhx <- list(
  down = rownames(significant_genes_stat_otf2_lhx[significant_genes_stat_otf2_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_otf2_lhx[significant_genes_stat_otf2_lhx$log2FoldChange > 0,])
)


diffreg_GOs_idvm_lhx <- getGOs(
  genelist = list_diffregs_forGO_idvm_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_idvm_pitx <- getGOs(
  genelist = list_diffregs_forGO_idvm_pitx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_lep_lhx <- getGOs(
  genelist = list_diffregs_forGO_lep_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_lep_pitx <- getGOs(
  genelist = list_diffregs_forGO_lep_pitx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_pgrn1_lhx <- getGOs(
  genelist = list_diffregs_forGO_pgrn1_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_pgrn1_pitx <- getGOs(
  genelist = list_diffregs_forGO_pgrn1_pitx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_tph_lhx <- getGOs(
  genelist = list_diffregs_forGO_tph_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_tph_pitx <- getGOs(
  genelist = list_diffregs_forGO_tph_pitx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_otf1_lhx <- getGOs(
  genelist = list_diffregs_forGO_otf1_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_otf2_lhx <- getGOs(
  genelist = list_diffregs_forGO_otf2_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)


pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/IDVM_LHX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_idvm_lhx$GOplot$down + ggtitle("IDVM_LHX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/IDVM_LHX_UP.pdf", width = 10, height = 10)
diffreg_GOs_idvm_lhx$GOplot$up + ggtitle("IDVM_LHX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/IDVM_PITX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_idvm_pitx$GOplot$down + ggtitle("IDVM_PITX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/IDVM_PITX_UP.pdf", width = 10, height = 10)
diffreg_GOs_idvm_pitx$GOplot$up + ggtitle("IDVM_PITX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/LEP_LHX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_lep_lhx$GOplot$down + ggtitle("LEP_LHX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/LEP_LHX_UP.pdf", width = 10, height = 10)
diffreg_GOs_lep_lhx$GOplot$up + ggtitle("LEP_LHX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/LEP_PITX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_lep_pitx$GOplot$down + ggtitle("LEP_PITX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/LEP_PITX_UP.pdf", width = 10, height = 10)
diffreg_GOs_lep_pitx$GOplot$up + ggtitle("LEP_PITX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/PNRG1_LHX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_pgrn1_lhx$GOplot$down + ggtitle("PNRG1_LHX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/PNRG1_LHX_UP.pdf", width = 10, height = 10)
diffreg_GOs_pgrn1_lhx$GOplot$up + ggtitle("PNRG1_LHX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/PNRG1_PITX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_pgrn1_pitx$GOplot$down + ggtitle("PNRG1_PITX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/PNRG1_PITX_UP.pdf", width = 10, height = 10)
diffreg_GOs_pgrn1_pitx$GOplot$up + ggtitle("PNRG1_PITX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/TPH_LHX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_tph_lhx$GOplot$down + ggtitle("TPH_LHX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/TPH_LHX_UP.pdf", width = 10, height = 10)
diffreg_GOs_tph_lhx$GOplot$up + ggtitle("TPH_LHX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/TPH_PITX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_tph_pitx$GOplot$down + ggtitle("TPH_PITX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/TPH_PITX_UP.pdf", width = 10, height = 10)
diffreg_GOs_tph_pitx$GOplot$up + ggtitle("TPH_PITX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/OTF1_LHX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_otf1_lhx$GOplot$down + ggtitle("OTF1_LHX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/OTF1_LHX_UP.pdf", width = 10, height = 10)
diffreg_GOs_otf1_lhx$GOplot$up + ggtitle("OTF1_LHX_UP")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/OTF2_LHX_DOWN.pdf", width = 10, height = 10)
diffreg_GOs_otf2_lhx$GOplot$down+ ggtitle("OTF2_LHX_DOWN")
dev.off()

pdf("/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/GOs/OTF2_LHX_UP.pdf", width = 10, height = 10)
diffreg_GOs_otf2_lhx$GOplot$up + ggtitle("OTF2_LHX_UP")
dev.off()



#Save Tables

head(rosetta)




significant_tph_pitx
significant_genes_stat_tph_pitx
significant_genes_tph_pitx <- rosetta[rosetta$gene %in% rownames(significant_tph_pitx), c(1, 6, 10, 11)]
head(significant_genes_tph_pitx)

table_tph_pitx <- merge(
  significant_genes_stat_tph_pitx,
  significant_genes_tph_pitx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_tph_pitx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_tph_pitx.tsv", sep = "\t", row.names = FALSE)


significant_tph_lhx
significant_genes_stat_tph_lhx
significant_genes_tph_lhx <- rosetta[rosetta$gene %in% rownames(significant_tph_lhx), c(1, 6, 10, 11)]
head(significant_tph_lhx)

table_tph_lhx <- merge(
  significant_genes_stat_tph_lhx,
  significant_genes_tph_lhx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_tph_lhx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_tph_lhx.tsv", sep = "\t", row.names = FALSE)


significant_otf2_lhx
significant_genes_stat_otf2_lhx
significant_genes_otf2_lhx <- rosetta[rosetta$gene %in% rownames(significant_otf2_lhx), c(1, 6, 10, 11)]
head(significant_otf2_lhx)

table_otf2_lhx <- merge(
  significant_genes_stat_otf2_lhx,
  significant_genes_otf2_lhx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_otf2_lhx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_otf2_lhx.tsv", sep = "\t", row.names = FALSE)


significant_pgrn1_lhx
significant_genes_stat_pgrn1_lhx
significant_genes_pgrn1_lhx <- rosetta[rosetta$gene %in% rownames(significant_pgrn1_lhx), c(1, 6, 10, 11)]
head(significant_pgrn1_lhx)

table_pgrn1_lhx <- merge(
  significant_genes_stat_pgrn1_lhx,
  significant_genes_pgrn1_lhx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_pgrn1_lhx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_pgrn1_lhx.tsv", sep = "\t", row.names = FALSE)


significant_pgrn1_pitx
significant_genes_stat_pgrn1_pitx
significant_genes_pgrn1_pitx <- rosetta[rosetta$gene %in% rownames(significant_pgrn1_pitx), c(1, 6, 10, 11)]
head(significant_genes_pgrn1_pitx)

table_pgrn1_pitx <- merge(
  significant_genes_stat_pgrn1_pitx,
  significant_genes_pgrn1_pitx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_pgrn1_pitx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_pgrn1_pitx.tsv", sep = "\t", row.names = FALSE)


significant_lep_lhx
significant_genes_stat_lep_lhx
significant_genes_lep_lhx <- rosetta[rosetta$gene %in% rownames(significant_lep_lhx), c(1, 6, 10, 11)]
head(significant_lep_lhx)

table_lep_lhx <- merge(
  significant_genes_stat_lep_lhx,
  significant_genes_lep_lhx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_lep_lhx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_lep_lhx.tsv", sep = "\t", row.names = FALSE)


significant_lep_pitx
significant_genes_stat_lep_pitx
significant_genes_lep_pitx <- rosetta[rosetta$gene %in% rownames(significant_lep_pitx), c(1, 6, 10, 11)]
head(significant_genes_lep_pitx)

table_lep_pitx <- merge(
  significant_genes_stat_lep_pitx,
  significant_genes_lep_pitx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_lep_pitx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_lep_pitx.tsv", sep = "\t", row.names = FALSE)



significant_idvm_lhx
significant_genes_stat_idvm_lhx
significant_genes_idvm_lhx <- rosetta[rosetta$gene %in% rownames(significant_idvm_lhx), c(1, 6, 10, 11)]
head(significant_idvm_lhx)

table_idvm_lhx <- merge(
  significant_genes_stat_idvm_lhx,
  significant_genes_idvm_lhx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_idvm_lhx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_idvm_lhx.tsv", sep = "\t", row.names = FALSE)


significant_idvm_pitx
significant_genes_stat_idvm_pitx
significant_genes_idvm_pitx <- rosetta[rosetta$gene %in% rownames(significant_idvm_pitx), c(1, 6, 10, 11)]
head(significant_genes_idvm_pitx)

table_idvm_pitx <- merge(
  significant_genes_stat_idvm_pitx,
  significant_genes_idvm_pitx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_idvm_pitx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/cl/table_idvm_pitx.tsv", sep = "\t", row.names = FALSE)

