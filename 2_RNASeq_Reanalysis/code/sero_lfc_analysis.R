# WHOLE PSEUDOBULK DGE ANALYSIS

# Setup
setwd("~/colabos/elena_sero/rnaseq_reanalysis")

dir <- "/mnt/sda/elena/R_analysis/sero_dge_202310/sero_DGE_per_cluser/"

source("/mnt/sda/alberto/colabos/sizes/code/R_code/r_functions/sourcefolder.R")
sourceFolder("/mnt/sda/alberto/colabos/sizes/code/R_code/r_functions/",recursive = TRUE)

f_col <- function(x){
  
  a = as.numeric(x[2])
  b = as.numeric(x[3])
  
  if(a < -1 & b < -1){
    y = ggplot2::alpha("lightblue",0.3)
  } else if(a > 1 & b > 1){
    y = ggplot2::alpha("salmon",0.3)
  } else{
    y = rgb(0,0,0,0.1)
  }
  
  return(y)
}

library(Matrix)
library(topGO)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggvenn)

# Load
load("outputs/rda/currie_et_al_reanalysis_lhx.rda")
load("outputs/rda/currie_et_al_reanalysis_pitx.rda")
load("outputs/rda/marz_et_al_reanalysis.rda")

# Load necessary data
sero_Idents <- read.delim2(paste0(dir,"06092023_smed_sero_identities.csv"), sep = "," , header = TRUE)
sero_genes <- read.delim2(paste0(dir,"06092023_smed_sero_genes.csv"), sep = "," , header = TRUE)
sero_leiden_col <- read.delim2(paste0(dir,"leiden_3_colors.csv"), sep = ",")
sero_Idents$color <- translate_ids(x=sero_Idents$leiden_3,dict = sero_leiden_col)
sero_Idents$whole <- "whole"


sero_ctypes <- unique(sero_Idents[,c("leiden_3","leiden_3_names","broad_names","color")])
sero_ctypes <- merge(
  sero_ctypes,
  as.data.frame(table(sero_Idents$leiden_3)),
  by.x = 1,
  by.y = 1,
  all.x = TRUE
)

rosetta <- read.delim2("/mnt/sda/alberto/projects/smed_rink_gene_annot/20231103_Smed_Rink_Simplified_Annotation_Table.tsv",header = TRUE,sep="\t")

# Load the matrix
sero_X <- readMM(file = paste0(dir,"matrix.mtx"))
colnames(sero_X) <-  read.table(paste0(dir,"barcodes.tsv"))[,1]
rownames(sero_X) <-  read.table(paste0(dir,"features.tsv"))[,1]
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

sero_pseudowhole_cond_reps <- pseudobulk_cond_rep(
  x = sero_X,
  identities = sero_Idents$whole,
  conditions = sero_Idents$sample,
  replicates = sero_Idents$group
)


# Clean sampletable

sero_sampletable <- sero_pseudowhole_cond_reps$sampletable

sero_sampletable_lhx = sero_sampletable[sero_sampletable$condition %in% c("gfp(RNAi)", "lhx1/5-1(RNAi)"), ] # Change based on your conditions
sero_sampletable_lhx <- clean_sampletable(sero_sampletable_lhx)

sero_sampletable_pitx = sero_sampletable[sero_sampletable$condition %in% c("gfp(RNAi)", "pitx(RNAi)"), ] # Change based on your conditions
sero_sampletable_pitx <- clean_sampletable(sero_sampletable_pitx)

sero_matrix_lhx <- 
  sero_pseudowhole_cond_reps$matrix[
    # rownames(sero_pseudowhole_cond_reps$matrix) %in% rosetta$gene[rosetta$gene_type == "hconf"],
    ,
    colnames(sero_pseudowhole_cond_reps$matrix) %in% sero_sampletable_lhx$sample
  ]

sero_matrix_pitx <- 
  sero_pseudowhole_cond_reps$matrix[
    # rownames(sero_pseudowhole_cond_reps$matrix) %in% rosetta$gene[rosetta$gene_type == "hconf"],
    ,
    colnames(sero_pseudowhole_cond_reps$matrix) %in% sero_sampletable_pitx$sample
  ]


# SINGLE CELL ANALYSIS OVER ALL CELL TYPES
sero_DGE_whole_lhx<-
  deseq_pseudobulk(
    count_matrix = sero_matrix_lhx,
    samples_info = sero_sampletable_lhx[,-2],
    celltype = "whole",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","lhx1/5-1(RNAi)","gfp(RNAi)"),
    plot_results = FALSE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )

sero_DGE_whole_pitx <-
  deseq_pseudobulk(
    count_matrix = sero_matrix_pitx,
    samples_info = sero_sampletable_pitx[,-2],
    celltype = "whole",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","pitx(RNAi)","gfp(RNAi)"),
    plot_results = FALSE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )


# LFC LHX
a = as.data.frame(res_lhx); a$gene <- as.character(rownames(a))
b = as.data.frame(sero_DGE_whole_lhx$res); b$gene <- as.character(rownames(b))
  
lhx_lfc <-
  merge(
    a[,c(7,2)],
    b[,c(7,2)],
    by=1,
    all = TRUE
  )
lhx_lfc <- lhx_lfc[complete.cases(lhx_lfc),]
colnames(lhx_lfc) <- c("gene","currie","emili")

plot(lhx_lfc$currie,lhx_lfc$emili)


# LFC LHX
a = as.data.frame(res_lhx); a$gene <- as.character(rownames(a))
b = as.data.frame(sero_DGE_whole_lhx$res); b$gene <- as.character(rownames(b))

lhx_lfc <-
  merge(
    a[,c(7,2)],
    b[,c(7,2)],
    by=1,
    all = TRUE
  )
lhx_lfc <- lhx_lfc[complete.cases(lhx_lfc),]
colnames(lhx_lfc) <- c("gene","currie","emili")

plot(
  lhx_lfc$currie,lhx_lfc$emili,
  main = "correlation LFCs (Lhx RNAi),\nCurrie et al 2013\nvs this study",
  xlab = "logFC in Currie et al., 2013",
  ylab = "logFC in this study"
  )




# LFC PITX
a = as.data.frame(res_pitx); a$gene <- as.character(rownames(a))
b = as.data.frame(sero_DGE_whole_pitx$res); b$gene <- as.character(rownames(b))

pitx_lfc <-
  merge(
    a[,c(7,2)],
    b[,c(7,2)],
    by=1,
    all = TRUE
  )
pitx_lfc <- pitx_lfc[complete.cases(pitx_lfc),]
colnames(pitx_lfc) <- c("gene","currie","emili")
plot(pitx_lfc$currie,pitx_lfc$emili)
plot(
  pitx_lfc$currie,pitx_lfc$emili,
  main = "correlation LFCs (Pitx RNAi),\nCurrie et al 2013\nvs this study",
  xlab = "logFC in Currie et al., 2013",
  ylab = "logFC in this study"
  )


# LFC PITX MARZ
a = as.data.frame(res); a$gene <- as.character(rownames(a))
b = as.data.frame(sero_DGE_whole_pitx$res); b$gene <- as.character(rownames(b))

pitx_lfc_marz <-
  merge(
    a[,c(7,2)],
    b[,c(7,2)],
    by=1,
    all = TRUE
  )
pitx_lfc_marz <- pitx_lfc_marz[complete.cases(pitx_lfc_marz),]
colnames(pitx_lfc_marz) <- c("gene","marz","emili")

# LFC PITX MARZ
a = as.data.frame(res); a$gene <- as.character(rownames(a))
b = as.data.frame(res_pitx); b$gene <- as.character(rownames(b))

pitx_currie_marz <-
  merge(
    a[,c(7,2)],
    b[,c(7,2)],
    by=1,
    all = TRUE
  )
pitx_currie_marz <- pitx_currie_marz[complete.cases(pitx_currie_marz),]
colnames(pitx_currie_marz) <- c("gene","marz","currie")

## SAVE PLOTS, TABLES

pdf("graphics/pitx_marx_currie_comparison.pdf")
plot(
  pitx_currie_marz$marz,pitx_currie_marz$currie,
  main = "correlation LFCs (Pitx RNAi),\nMarz et al 2013\nvs Currie et al 2013",
  xlab = "logFC in Marz et al., 2013",
  ylab = "logFC in Cuerie et al., 2013",
  ylim = c(-5,5),
  pch = 21,
  col = NA,
  bg = apply(pitx_currie_marz,1,f_col)
)
abline(v=c(-1,1), lty = 2)
abline(h=c(-1,1), lty = 2)
dev.off()

pdf("graphics/pitx_marx_lfc_comparison.pdf")
plot(
  pitx_lfc_marz$marz,pitx_lfc_marz$emili,
  main = "correlation LFCs (Pitx RNAi),\nMarz et al 2013\nvs this study",
  xlab = "logFC in Marz et al., 2013",
  ylab = "logFC in this study",
  ylim = c(-5,5),
  pch = 21,
  col = NA,
  bg = apply(pitx_lfc_marz,1,f_col)
)
abline(v=c(-1,1), lty = 2)
abline(h=c(-1,1), lty = 2)
dev.off()

pdf("graphics/pitx_currie_lfc_comparison.pdf")
plot(
  pitx_lfc$currie,pitx_lfc$emili,
  main = "correlation LFCs (Pitx RNAi),\nCurrie et al 2013\nvs this study",
  xlab = "logFC in Currie et al., 2013",
  ylab = "logFC in this study",
  ylim = c(-5,5),
  xlim = c(-5,5),
  pch = 21,
  col = NA,
  bg = apply(pitx_lfc,1,f_col)
)
abline(v=c(-1,1), lty = 2)
abline(h=c(-1,1), lty = 2)
dev.off()

pdf("graphics/lhx_lfc_comparison.pdf")
plot(
  lhx_lfc$currie,lhx_lfc$emili,
  main = "correlation LFCs (Lhx RNAi),\nCurrie et al 2013\nvs this study",
  xlab = "logFC in Currie et al., 2013",
  ylab = "logFC in this study",
  ylim = c(-5,5),
  xlim = c(-5,5),
  pch = 21,
  col = NA,
  bg = apply(lhx_lfc,1,f_col)
)
abline(v=c(-1,1), lty = 2)
abline(h=c(-1,1), lty = 2)
dev.off()

negFC_marz_emili <- 
  rosetta[
    rosetta$gene %in%
      pitx_lfc_marz$gene[
        which(pitx_lfc_marz$marz < -1 &
          pitx_lfc_marz$emili < -1)
      ],
  ]

posFC_marz_emili <-
  rosetta[
    rosetta$gene %in%
      pitx_lfc_marz$gene[
        which(pitx_lfc_marz$marz > 1 &
                pitx_lfc_marz$emili > 1)
      ],
  ]

write.table(
  negFC_marz_emili,
  "outputs/neg_FC_marz_emili.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  posFC_marz_emili,
  "outputs/pos_FC_marz_emili.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

negFC_currie_emili_lhx <- 
  rosetta[
    rosetta$gene %in%
      lhx_lfc$gene[
        which(lhx_lfc$currie < -1 &
                lhx_lfc$emili < -1)
      ],
  ]

posFC_currie_emili_lhx <-
  rosetta[
    rosetta$gene %in%
      lhx_lfc$gene[
        which(lhx_lfc$currie > 1 &
                lhx_lfc$emili > 1)
      ],
  ]

write.table(
  negFC_currie_emili_lhx,
  "outputs/negFC_currie_emili_lhx.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  posFC_currie_emili_lhx,
  "outputs/posFC_currie_emili_lhx.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

negFC_currie_emili_pitx <- 
  rosetta[
    rosetta$gene %in%
      pitx_lfc$gene[
        which(pitx_lfc$currie < -1 &
                pitx_lfc$emili < -1)
      ],
  ]

posFC_currie_emili_pitx <-
  rosetta[
    rosetta$gene %in%
      pitx_lfc$gene[
        which(pitx_lfc$currie > 1 &
                pitx_lfc$emili > 1)
      ],
  ]

write.table(
  negFC_currie_emili_pitx,
  "outputs/negFC_currie_emili_pitx.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  posFC_currie_emili_pitx,
  "outputs/posFC_currie_emili_pitx.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

out <- 
  merge(
    data.frame(
      gene = lhx_lfc$gene,
      Currie_et_al_Lhx_LFC = lhx_lfc$currie,
      This_study_Lhx_LFC = lhx_lfc$emili
    ),
    data.frame(
      gene = pitx_lfc$gene,
      Currie_et_al_Pitx_LFC = pitx_lfc$currie,
      This_study_Litx_LFC = pitx_lfc$emili
    ),
    by = 1,
    all = TRUE
  )

write.table(
  out,
  file = "outputs/logfoldchanges_currie_and_this_study.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  dec = ".",
  quote = FALSE
)
