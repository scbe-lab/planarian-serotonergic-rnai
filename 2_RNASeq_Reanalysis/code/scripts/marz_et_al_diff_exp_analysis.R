library(DESeq2)
library(tximport)
library(ComplexHeatmap)
library(circlize)
library(colorspace)

dir <- "~/colabos/elena_sero/rnaseq_reanalysis"
setwd(dir)

samplesdir <-
  "outputs/RNA/kallisto_output/raw_output/"

sampleTable <-
  read.table(
    "outputs/RNA/kallisto_output/sampletable.tsv",
    header = F
  )

colnames(sampleTable) <-
  c("name", "condition", "batch", "dir")

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

res <- results(object = dds, contrast = c("condition","pitx_RNAi","GFP_RNAi"))
r <- as.data.frame(res)
r <-r [complete.cases(r),]

# Retrieve diff genes
filter_by = "padj"
p_threshold = 0.1
filter_by_column <- which(colnames(r) == filter_by)
signif <- rownames(r)[r[,filter_by_column] < p_threshold]
signif <- signif[complete.cases(signif)]


v <-
  r %>%
  ggplot(
    aes(x = log2FoldChange, y = -log(!! sym(filter_by)), color = factor(!! sym(filter_by) < p_threshold))
  ) +
  geom_point() +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red3")) +
  guides(color = "none") + 
  theme_minimal()
v

# Create transformed data matrix for heatmap
print("heatmap qnorm")
m <- counts(dds,normalize = TRUE)
vsd <-  varianceStabilizingTransformation(dds, blind = FALSE)
m_hm <- assay(vsd[rownames(vsd) %in% signif,])
colnames(m_hm) <- gsub("kallisto_out_","",colnames(m_hm))
# m_hm <- m[rownames(m) %in% signif,]
# m_hm <-  as.matrix(t(scale(t(quantnorm(m)[rownames(m) %in% signif,]))))

# Create heatmap annotations
ha_sampletable <- HeatmapAnnotation(
  df = sampleTable[,-c(1,4)], 
  col = dynamic_colors_annotation(df = sampleTable[,-c(1,4)], rand.seed = 2)
)

# Create heatmap
hm <- Heatmap(
  name = "expression",
  t(scale(t(m_hm))),
  clustering_method_rows = "ward.D2",
  bottom_annotation = ha_sampletable,
  column_title = "pitx RNAi (Marz et al., 2013)",
  show_row_names = FALSE,
  col = diverging_hcl(10,"Blue-Red 2")
)
hm

# Create pca
pcaData <- plotPCA(vsd, intgroup=c("condition", "batch"), returnData=TRUE)
pcaData$batch <- factor(pcaData$batch, levels = c("1","2"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
pc <- 
  ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
pc


signif_data <- rosetta[rosetta$gene %in% signif,]

save(
  sampleTable, # table with experiment/samples/conditions info
  dds, # deseq2 object, already al calculations done
  res, # results table, unfiltered
  r, # results table, filtered significant p.adj < 0.1
  v, # volcano plot
  hm, # heatmap
  pc, # PCA
  signif, # list of DEG genes
  signif_data, # rosetta of DEG genes
  file = "outputs/rda/marz_et_al_reanalysis.rda"
)