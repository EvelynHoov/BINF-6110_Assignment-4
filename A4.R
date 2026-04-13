library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(Matrix)


# normalization
head(seurat_ass4@meta.data)
seurat_ass4@assays$RNA@data[1:5, 1:5]
dim(seurat_ass4@assays$RNA@scale.data)
"pca" %in% names(seurat_ass4@reductions)
table(seurat_ass4$seurat_clusters)
"umap" %in% names(seurat_ass4@reductions)

seurat_ass4 <- FindVariableFeatures(seurat_ass4, selection.method = "vst", nfeatures = 2000)


seurat_ass4 <- ScaleData(seurat_ass4)


seurat_ass4 <- RunPCA(seurat_ass4)


seurat_ass4 <- FindNeighbors(seurat_ass4, dims = 1:30)


seurat_ass4 <- FindClusters(seurat_ass4, resolution = 0.30)


seurat_ass4 <- RunUMAP(seurat_ass4, dims = 1:30)
length(VariableFeatures(seurat_ass4))

# clustering
# Clusters on UMAP
DimPlot(seurat_ass4, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP of clusters")
# By sample group (replace 'sample_group' with e.g. 'time' or 'disease__ontology_label')
DimPlot(seurat_ass4, reduction = "umap", group.by = "organ_custom") +
  ggtitle("UMAP by sample group")

DefaultAssay(seurat_ass4) <- "RNA"
cl3 <- subset(seurat_ass4, idents = "3")

meta <- cl3@meta.data[, c("organ_custom", "mouse_id")]
meta$sample_id <- paste(meta$mouse_id, meta$organ_custom, sep = "_")

counts <- cl3@assays$RNA@counts

# aggregate counts by sample_id
pb_counts <- rowsum(t(counts), group = meta$sample_id)
pb_counts <- t(pb_counts)


sample_info <- data.frame(sample_id = colnames(pb_counts)) %>%
  tidyr::extract(
    sample_id,
    into = c("mouse_id", "organ_custom"),
    regex = "^(.*)_(LNG|OM|RM)$"
  )

# collapse OM + RM into OM_RM
sample_info$organ_group <- ifelse(sample_info$organ_custom == "LNG", "LNG", "OM_RM")
sample_info$organ_group <- factor(sample_info$organ_group, levels = c("OM_RM", "LNG"))

rownames(sample_info) <- sample_info$sample_id

dds <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData   = sample_info,
  design    = ~ organ_group
)

dds <- DESeq(dds)

res_LNG_vs_OMRM <- results(dds, contrast = c("organ_group", "LNG", "OM_RM"))
res_LNG_vs_OMRM <- as.data.frame(res_LNG_vs_OMRM)
res_LNG_vs_OMRM$gene <- rownames(res_LNG_vs_OMRM)

FeaturePlot(
  seurat_ass4,
  features = c("Obp1a", "Ifit1", "Fosl2", "Gnb1"),
  cols = c("lightgrey", "red"),
  reduction = "umap",
  ncol = 2
)

ggplot(res_LNG_vs_OMRM, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(
    title = "Cluster 3: LNG vs OM+RM (DESeq2)",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  )


sig <- res_LNG_vs_OMRM %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.5)

ego <- enrichGO(
  gene     = sig$gene,
  OrgDb    = org.Mm.eg.db,
  keyType  = "SYMBOL",
  ont      = "BP",
  readable = TRUE
)

dotplot(ego, showCategory = 20) +
  ggtitle("Cluster 3 ORA: LNG vs OM+RM")


gene_list <- res_LNG_vs_OMRM$log2FoldChange
names(gene_list) <- res_LNG_vs_OMRM$gene
gene_list <- sort(gene_list, decreasing = TRUE)

gsea <- gseGO(
  geneList = gene_list,
  OrgDb    = org.Mm.eg.db,
  keyType  = "SYMBOL",
  ont      = "BP"
)

dotplot(gsea, showCategory = 20) +
  ggtitle("Cluster 3 GSEA: LNG vs OM+RM")

gseaplot2(gsea, geneSetID = 1)