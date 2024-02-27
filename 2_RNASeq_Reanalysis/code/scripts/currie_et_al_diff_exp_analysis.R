library(DESeq2)
library(tximport)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(dplyr)
library(ggplot2)

source("~/colabos/sizes/code/R_code/r_functions/sourcefolder.R")
sourceFolder("~/colabos/sizes/code/R_code/r_functions/")

rosetta <- read.delim2("~/projects/smed_rink_gene_annot/20231103_Smed_Rink_Simplified_Annotation_Table.tsv")

dir <- "~/colabos/elena_sero/rnaseq_reanalysis"
setwd(dir)

samplesdir <-
  "outputs/RNA/output_currie/kallisto_output/raw_output/"

sampleTable <-
  read.table(
    "outputs/RNA/output_currie/sampletable.tsv",
    sep = "\t",
    header = TRUE
  )[,c("Run","condition","replicate","path")]

colnames(sampleTable) <-
  c("name", "condition", "batch", "dir")

sampleTable$condition <- factor(sampleTable$condition, levels = c("control(RNAi)", "lhx1/5-1(RNAi)", "pitx(RNAi)"))
sampleTable$batch <- factor(sampleTable$batch, levels = c("1","2","3"))
sampleTable

files <-
  setNames(
    object = list.files(
      samplesdir,
      "h5",
      recursive = T,
      full.names = T
    ),
    list.dirs(
      samplesdir,
      full.names = F
    )[-1]
  )

txi.kallisto <-
  tximport(
    files,
    type = "kallisto",
    txOut = TRUE
  )

dds <-
  DESeqDataSetFromTximport(
    txi.kallisto,
    colData = sampleTable,
    design = ~ batch + condition
  )

# convert transcript IDs to gene IDs
rownames(dds) <- translate_ids(x = rownames(dds), dict = rosetta[,c(2,1)])

dds <- DESeq(dds)

# Create pca
vsd <-  varianceStabilizingTransformation(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "batch"), returnData=TRUE)
pcaData$batch <- factor(pcaData$batch, levels = c("1","2","3"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
pc <- 
  ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
pc


# pick by which to filter
filter_by = "padj"
p_threshold = 0.05

## PITX RNAi
contrast_pitx <- c("pitx(RNAi)",  "control(RNAi)")
res_pitx <- results(object = dds, contrast =  c("condition", contrast_pitx))
r_pitx <- as.data.frame(res_pitx)
r_pitx <-r_pitx [complete.cases(r_pitx),]

# Retrieve diff genes
filter_by_column <- which(colnames(r_pitx) == filter_by)
signif_pitx <- rownames(r_pitx)[r_pitx[,filter_by_column] < p_threshold]
signif_pitx <- signif_pitx[complete.cases(signif_pitx)]


v_pitx <-
  r_pitx %>%
  ggplot(
    aes(x = log2FoldChange, y = -log(!! sym(filter_by)), color = factor(!! sym(filter_by) < p_threshold))
  ) +
  geom_point() +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red3")) +
  guides(color = "none") + 
  theme_minimal()
v_pitx

# Create transformed data matrix for heatmap
m_hm_pitx <- 
  assay(vsd[
    rownames(vsd) %in% signif_pitx,
    gsub("kallisto_out_", "", colnames(vsd)) %in% sampleTable$name[sampleTable$condition %in% contrast_pitx]
    ])
colnames(m_hm_pitx) <- gsub("kallisto_out_","",colnames(m_hm_pitx))

# Create heatmap annotations
ha_sampletable <- HeatmapAnnotation(
  df = sampleTable[sampleTable$condition %in% contrast_pitx,-c(1,4)], 
  col = dynamic_colors_annotation(df = sampleTable[sampleTable$condition %in% contrast_pitx,-c(1,4)], rand.seed = 2)
)

# Create heatmap
hm_pitx <- Heatmap(
  name = "expression",
  t(scale(t(m_hm_pitx))),
  clustering_method_rows = "ward.D2",
  bottom_annotation = ha_sampletable,
  column_title = "pitx RNAi (Currie et al., 2013)",
  show_row_names = FALSE,
  col = diverging_hcl(10,"Blue-Red 2")
)
hm_pitx

signif_data_pitx <- rosetta[rosetta$gene %in% signif_pitx,]

## LHX RNAi
contrast_lhx <- c("lhx1/5-1(RNAi)",  "control(RNAi)")
res_lhx <- results(object = dds, contrast =  c("condition", contrast_lhx))
r_lhx <- as.data.frame(res_lhx)
r_lhx <-r_lhx [complete.cases(r_lhx),]

# Retrieve diff genes
filter_by_column <- which(colnames(r_lhx) == filter_by)
signif_lhx <- rownames(r_lhx)[r_lhx[,filter_by_column] < p_threshold]
signif_lhx <- signif_lhx[complete.cases(signif_lhx)]


v_lhx <-
  r_lhx %>%
  ggplot(
    aes(x = log2FoldChange, y = -log(!! sym(filter_by)), color = factor(!! sym(filter_by) < p_threshold))
  ) +
  geom_point() +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red3")) +
  guides(color = "none") + 
  theme_minimal()
v_lhx

# Create transformed data matrix for heatmap
m_hm_lhx <- 
  assay(vsd[
    rownames(vsd) %in% signif_lhx,
    gsub("kallisto_out_", "", colnames(vsd)) %in% sampleTable$name[sampleTable$condition %in% contrast_lhx]
  ])
colnames(m_hm_lhx) <- gsub("kallisto_out_","",colnames(m_hm_lhx))

# Create heatmap annotations
ha_sampletable <- HeatmapAnnotation(
  df = sampleTable[sampleTable$condition %in% contrast_lhx,-c(1,4)], 
  col = dynamic_colors_annotation(df = sampleTable[sampleTable$condition %in% contrast_lhx,-c(1,4)], rand.seed = 2)
)

# Create heatmap
hm_lhx <- Heatmap(
  name = "expression",
  t(scale(t(m_hm_lhx))),
  clustering_method_rows = "ward.D2",
  bottom_annotation = ha_sampletable,
  column_title = "lhx1/5-1 RNAi (Currie et al., 2013)",
  show_row_names = FALSE,
  col = diverging_hcl(10,"Blue-Red 2")
)
hm_lhx

signif_data_lhx <- rosetta[rosetta$gene %in% signif_lhx,]

# Save data
save(
  sampleTable, # table with experiment/samples/conditions info
  dds, # deseq2 object, already al calculations done
  pc, # PCA
  res_pitx, # results table, unfiltered
  r_pitx, # results table, filtered significant p.adj < 0.1
  v_pitx, # volcano plot
  hm_pitx, # heatmap
  signif_pitx, # list of DEG genes
  signif_data_pitx, # rosetta of DEG genes
  file = "outputs/rda/currie_et_al_reanalysis_pitx.rda"
)

save(
  sampleTable, # table with experiment/samples/conditions info
  dds, # deseq2 object, already al calculations done
  pc, # PCA
  res_lhx, # results table, unfiltered
  r_lhx, # results table, filtered significant p.adj < 0.1
  v_lhx, # volcano plot
  hm_lhx, # heatmap
  signif_lhx, # list of DEG genes
  signif_data_lhx, # rosetta of DEG genes
  file = "outputs/rda/currie_et_al_reanalysis_lhx.rda"
)