# SERO DGE ANALYSIS

# Setup
dir <- "/mnt/sda/elena/R_analysis/sero_dge_202310/sero_DGE_per_broadgroup/"
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
sero_broad_col <- read.delim2("broad_colors.csv", sep = ",")
sero_Idents$color <- translate_ids(x=sero_Idents$broad_names,dict = sero_broad_col)

sero_ctypes <- unique(sero_Idents[,c("broad_names","color")])
sero_ctypes <- merge(
  sero_ctypes,
  as.data.frame(table(sero_Idents$broad_names)),
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
  identities = sero_Idents$broad_names,
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
    by.y = "broad_names",
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
    by.y = "broad_names",
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

pseudobulk_matrix <- pseudobulk(x = as.matrix(sero_X), identities = sero_Idents$broad_names)

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

sero_DGE_m_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "muscle",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
m_lhx <- sero_DGE_m_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "muscle DEGs lhx1/5-1(RNAi)")
m_lhx


sero_DGE_m_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "muscle",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
m_pitx <- sero_DGE_m_pitx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "muscle DEGs pitx(RNAi)")
m_pitx


sero_DGE_e_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "epidermal",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
e_lhx <- sero_DGE_e_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "epidermal DEGs lhx1/5-1(RNAi)")
e_lhx


sero_DGE_e_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "epidermal",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
e_pitx <- sero_DGE_e_pitx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "epidermal DEGs pitx(RNAi)")
e_pitx


sero_DGE_p_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "parenchymal",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
p_lhx <- sero_DGE_p_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "parenchymal DEGs lhx1/5-1(RNAi)")
p_lhx


sero_DGE_p_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "parenchymal",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
p_pitx <- sero_DGE_p_pitx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "parenchymal DEGs pitx(RNAi)")
p_pitx


sero_DGE_n_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "neurons",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
n_pitx <- sero_DGE_n_pitx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "gray", fill=NA, linewidth = 0.75))+
  ggtitle(label = "neurons DEGs pitx(RNAi)")
n_pitx


sero_DGE_n_lhx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "neurons",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
n_lhx <- sero_DGE_n_lhx$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "neurons DEGs lhx1/5-1(RNAi)")
n_lhx


sero_DGE_m_lhx$heatmap
sero_DGE_m_pitx$heatmap

sero_DGE_e_lhx$heatmap
sero_DGE_e_pitx$heatmap


sero_DGE_p_lhx$heatmap
sero_DGE_p_pitx$heatmap

sero_DGE_n_pitx$heatmap
sero_DGE_n_lhx$heatmap



# GO enrichment analysis

source("/mnt/sda/alberto/colabos/sizes/code/R_code/r_functions/sourcefolder.R")
sourceFolder("/mnt/sda/alberto/colabos/sizes/code/R_code/r_functions/",recursive = TRUE)
source("/mnt/sda/alberto/projects/smed_cisreg/code/r_code/functions/topGO_wrapper.R")

smed_id_GO <- readMappings("/mnt/sda/alberto/projects/smed_cisreg/outputs/gene_annotation/smed_GOs.tsv")

sero_DGE_m_lhx$res
sero_DGE_m_pitx$res
sero_DGE_e_lhx$res
sero_DGE_e_pitx$res
sero_DGE_p_lhx$res
sero_DGE_p_pitx$res
sero_DGE_n_pitx$res
sero_DGE_n_lhx$res


significant_m_lhx <- subset(sero_DGE_m_lhx$res, pvalue < 0.05)
significant_m_pitx <- subset(sero_DGE_m_pitx$res, pvalue < 0.05)
significant_e_lhx <- subset(sero_DGE_e_lhx$res, pvalue < 0.05)
significant_e_pitx <- subset(sero_DGE_e_pitx$res, pvalue < 0.05)
significant_p_lhx <- subset(sero_DGE_p_lhx$res, pvalue < 0.05)
significant_p_pitx <- subset(sero_DGE_p_pitx$res, pvalue < 0.05)
significant_n_lhx <- subset(sero_DGE_n_lhx$res, pvalue < 0.05)
significant_n_pitx <- subset(sero_DGE_n_pitx$res, pvalue < 0.05)

significant_genes_stat_m_lhx <- data.frame(significant_m_lhx)
significant_genes_stat_m_pitx <- data.frame(significant_m_pitx)
significant_genes_stat_e_lhx <- data.frame(significant_e_lhx)
significant_genes_stat_e_pitx <- data.frame(significant_e_pitx)
significant_genes_stat_p_lhx <- data.frame(significant_p_lhx)
significant_genes_stat_p_pitx <- data.frame(significant_p_pitx)
significant_genes_stat_n_lhx <- data.frame(significant_n_lhx)
significant_genes_stat_n_pitx <- data.frame(significant_n_pitx)



list_diffregs_forGO_m_lhx <- list(
  down = rownames(significant_genes_stat_m_lhx[significant_genes_stat_m_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_m_lhx[significant_genes_stat_m_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_m_pitx <- list(
  down = rownames(significant_genes_stat_m_pitx[significant_genes_stat_m_pitx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_m_pitx[significant_genes_stat_m_pitx$log2FoldChange > 0,])
)
list_diffregs_forGO_e_lhx <- list(
  down = rownames(significant_genes_stat_e_lhx[significant_genes_stat_e_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_e_lhx[significant_genes_stat_e_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_e_pitx <- list(
  down = rownames(significant_genes_stat_e_pitx[significant_genes_stat_e_pitx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_e_pitx[significant_genes_stat_e_pitx$log2FoldChange > 0,])
)
list_diffregs_forGO_p_lhx <- list(
  down = rownames(significant_genes_stat_p_lhx[significant_genes_stat_p_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_p_lhx[significant_genes_stat_p_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_p_pitx <- list(
  down = rownames(significant_genes_stat_p_pitx[significant_genes_stat_p_pitx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_p_pitx[significant_genes_stat_p_pitx$log2FoldChange > 0,])
)
list_diffregs_forGO_n_lhx <- list(
  down = rownames(significant_genes_stat_n_lhx[significant_genes_stat_n_lhx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_n_lhx[significant_genes_stat_n_lhx$log2FoldChange > 0,])
)
list_diffregs_forGO_n_pitx <- list(
  down = rownames(significant_genes_stat_n_pitx[significant_genes_stat_n_pitx$log2FoldChange < 0,]),
  up = rownames(significant_genes_stat_n_pitx[significant_genes_stat_n_pitx$log2FoldChange > 0,])
)



diffreg_GOs_m_lhx <- getGOs(
  genelist = list_diffregs_forGO_m_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_m_pitx <- getGOs(
  genelist = list_diffregs_forGO_m_pitx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_e_lhx <- getGOs(
  genelist = list_diffregs_forGO_e_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_e_pitx <- getGOs(
  genelist = list_diffregs_forGO_e_pitx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_p_lhx <- getGOs(
  genelist = list_diffregs_forGO_p_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_p_pitx <- getGOs(
  genelist = list_diffregs_forGO_p_pitx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_n_lhx <- getGOs(
  genelist = list_diffregs_forGO_n_lhx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)
diffreg_GOs_n_pitx <- getGOs(
  genelist = list_diffregs_forGO_n_pitx,
  gene_universe= rownames(sero_X),
  gene2GO = smed_id_GO
)



diffreg_GOs_m_lhx$GOplot$down
diffreg_GOs_m_lhx$GOplot$up

diffreg_GOs_m_pitx$GOplot$down
diffreg_GOs_m_pitx$GOplot$up

diffreg_GOs_e_lhx$GOplot$down
diffreg_GOs_e_lhx$GOplot$up

diffreg_GOs_e_pitx$GOplot$down
diffreg_GOs_e_pitx$GOplot$up


diffreg_GOs_p_lhx$GOplot$down
diffreg_GOs_p_lhx$GOplot$up

diffreg_GOs_p_pitx$GOplot$down
diffreg_GOs_p_pitx$GOplot$up


diffreg_GOs_n_lhx$GOplot$down
diffreg_GOs_n_lhx$GOplot$up

diffreg_GOs_n_pitx$GOplot$down
diffreg_GOs_n_pitx$GOplot$up


save(sero_DGE_all_lhx, sero_DGE_all_pitx, file = "broad_sero_dge.Rda")


#Save Tables

head(rosetta)

significant_e_lhx
significant_genes_stat_e_lhx
significant_genes_e_lhx <- rosetta[rosetta$gene %in% rownames(significant_e_lhx), c(1, 6, 10, 11)]
head(significant_genes_e_lhx)

table_e_lhx <- merge(
  significant_genes_stat_e_lhx,
  significant_genes_e_lhx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
  )

write.table(table_e_lhx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/broad/table_e_lhx.tsv", sep = "\t", row.names = FALSE)


significant_n_pitx
significant_genes_stat_n_pitx
significant_genes_n_pitx <- rosetta[rosetta$gene %in% rownames(significant_n_pitx), c(1, 6, 10, 11)]
head(significant_genes_n_pitx)

table_n_pitx <- merge(
  significant_genes_stat_n_pitx,
  significant_genes_n_pitx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_n_pitx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/broad/table_n_pitx.tsv", sep = "\t", row.names = FALSE)


significant_n_lhx
significant_genes_stat_n_lhx
significant_genes_n_lhx <- rosetta[rosetta$gene %in% rownames(significant_n_lhx), c(1, 6, 10, 11)]
head(significant_genes_n_lhx)

table_n_lhx <- merge(
  significant_genes_stat_n_lhx,
  significant_genes_n_lhx,
  by.x = 0,
  by.y = 1,
  all.x = TRUE
)

write.table(table_n_lhx, "/mnt/sda/elena/sero_manuscript/DGE/graphic/broad/table_n_lhx.tsv", sep = "\t", row.names = FALSE)


