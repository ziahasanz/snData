# Zia - April 13th july, 2022

setwd("singlecell")
library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)

data <- Read10X("../singlecell/raw_feature_bc_matrix")
s1f <- CreateSeuratObject(counts = data,project = 'hippo', min.cells =1, min.features=20)

library(tidyverse)
library(Matrix)
library(cowplot)
theme_set(theme_cowplot())

s1f <- PercentageFeatureSet(s1f, pattern = "^Mt-", col.name = "percent.Mt")
VlnPlot(s1f, features =c('nFeature_RNA','nCount_RNA','percent.Mt'), pt.size =0.1)


s1f.sub <- subset(s1f, subset = nFeature_RNA > 40 & nCount_RNA > 400 & percent.Mt < 5)
VlnPlot(s1f, features = c("nFeature_RNA","nCount_RNA","percent.Mt"),ncol = 3)


s1f.sub.norm <- NormalizeData(s1f.sub, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)

#Identification of Highly Variable Features
s1f.sub.norm.vf <- FindVariableFeatures(s1f.sub.norm, selection.method = "vst", nfeatures = 2000,verbose = FALSE)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(s1f.sub.norm.vf), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(s1f.sub.norm.vf) + theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0,ynudge = 0)
plot1 + plot2

#scaling the data
all.genes <- rownames(s1f.sub.norm.vf)
s1f.sub.norm.vf.scale<- ScaleData(s1f.sub.norm.vf, features = all.genes, verbose = FALSE)
s1f.sub.norm.vf <- FindVariableFeatures(object = s1f.sub.norm.vf)

#if shows error: cannot allocate vector of size 4.5 Gb. Run the below command 
s1f.sub.norm.vf.scale <- ScaleData(s1f.sub.norm.vf)
s1f.sub.norm.vf.scale.PCA <- RunPCA(s1f.sub.norm.vf.scale)

# Examine and visualize PCA results a few different ways
print(s1f.sub.norm.vf.scale.PCA[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(s1f.sub.norm.vf.scale.PCA, dims = 1:2, reduction = "pca")
DimPlot(s1f.sub.norm.vf.scale.PCA, reduction = "pca")
DimHeatmap(s1f.sub.norm.vf.scale.PCA, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(s1f.sub.norm.vf.scale.PCA, dims = 1:15, cells = 500, balanced = TRUE)


#Determine the âdimensionalityâ of the dataset
#if you use below two commands for theâdimensionalityâ, it will take longtime to analyze the data (for instance: 4h for my data). I skip this step and use only the Elbow plot.
s1f.sub.norm.vf.scale.PCA.Dimenstion <- JackStraw(s1f.sub.norm.vf.scale.PCA, num.replicate = 100)
s1f.sub.norm.vf.scale.PCA.Dimenstion <- ScoreJackStraw(s1f.sub.norm.vf.scale.PCA.Dimenstion, dims = 1:7)
JackStrawPlot(s1f.sub.norm.vf.scale.PCA, dims = 1:15)
ElbowPlot(s1f.sub.norm.vf.scale.PCA)

#Cluster the cells
s1f.sub.norm.vf.scale.PCA.cluster <- FindNeighbors(s1f.sub.norm.vf.scale.PCA, dims = 1:15)
s1f.sub.norm.vf.scale.PCA.cluster <- FindClusters(s1f.sub.norm.vf.scale.PCA.cluster, resolution = 0.5)

#Adjust Resolution(0.3,0.5,0.7,1)  the higher the no., higher the clusters
s1f.sub.norm.vf.scale.PCA.cluster <- FindClusters(s1f.sub.norm.vf.scale.PCA.cluster, resolution = c(0.3,0.5,0.7,1))
View(s1f.sub.norm.vf.scale.PCA.cluster@meta.data)
DimPlot(s1f.sub.norm.vf.scale.PCA.cluster, group.by = "RNA_snn_res.1", label = TRUE)

#If you use resolution command then you have to set the identity of the clusters.Ser the "RNA_snn_res. which is suitable for the resolution 0.1-1.
Idents(s1f.sub.norm.vf.scale.PCA.cluster)
s1f.sub.norm.vf.scale.PCA.cluster <- "RNA_snn_res.0.5"

# Look at cluster IDs of the first 5 cells
head(Idents(s1f.sub.norm.vf.scale.PCA.cluster), 5)
s1f.sub.norm.vf.scale.PCA.cluster.umap <- RunUMAP(s1f.sub.norm.vf.scale.PCA.cluster, dims = 1:15)
DimPlot(s1f.sub.norm.vf.scale.PCA.cluster.umap, reduction = "umap")


#You can save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.
saveRDS(s1f.sub.norm.vf.scale.PCA.cluster.umap, file = "../singlecell/s1f_data_analysis.rds")
#tSNE Plot
s1f.sub.norm.vf.scale.PCA.cluster.tSNE <- RunTSNE(s1f.sub.norm.vf.scale.PCA.cluster, dims = 1:15, verbose = FALSE)
DimPlot(s1f.sub.norm.vf.scale.PCA.cluster.tSNE, reduction = "tsne")
DimPlot(s1f.sub.norm.vf.scale.PCA.cluster.tSNE, reduction = "tsne", label = TRUE)

#Finding differentially expressed features (cluster biomarkers)
#find all markers of cluster 0 (for individual cluster's marker)
cluster0.markers <- FindMarkers(object = s1f.sub.norm.vf.scale.PCA.cluster.umap, ident.1 = 1, min.pct = 0.25)
head(x = cluster0.markers, n = 10)
cluster0.markers <- FindMarkers(s1f.sub.norm.vf.scale.PCA.cluster.umap, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(s1f.sub.norm.vf.scale.PCA.cluster.umap, features = c("Zbtb20", "Map1b", "Erc2", "Olfm17", "Epha7", "Ppfia2", "Ppp3ca", "Rfx3",  "Adcy1"))
#this works
VlnPlot(s1f.sub.norm.vf.scale.PCA.cluster.umap, features = c("Gad1", "Gad2", "mog", "Slc17a7", "Pdgfra", "P2ry12", "Snap25", "Neurod6", "Megf11","PcdH15","Pdgfra", "Cspg4", "Slc1a2", "Slc1a3", "Glul", "S100b"))
VlnPlot(s1f.sub.norm.vf.scale.PCA.cluster.umap, features = c("Rbfox3", "Slc17a7", "Neurod6", "Snap25", "Gad1", "Gad2", "Megf11", "PcdH15","Pdgfra", "Cspg4", "P2ry12" , "Cx3cr1" , "Sall1" , "Csf1r", "Slc1a2", "Slc1a3", "Glul", "S100b", "Mog", "Mbp", "Fa2h", "Ugt8"))

# find markers for every cluster compared to all remaining cells, report only the positive ones
s1f.sub.norm.vf.scale.PCA.cluster.umap.markers <- FindAllMarkers(object = s1f.sub.norm.vf.scale.PCA.cluster.umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- s1f.sub.norm.vf.scale.PCA.cluster.umap.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

png(filename = "cluster_umap_heatmap.png", width = 1000, height = 1600, pointsize = 12, family = "sans", bg = "white")
DoHeatmap(s1f.sub.norm.vf.scale.PCA.cluster.umap, features = top10$gene) + NoLegend() + scale_fill_viridis(option = "B") + theme(text = element_text(size = 20))
dev.off()

## rank and save the marker genes

s1f.sub.norm.vf.scale.PCA.cluster.umap.markers$rank <-(s1f.sub.norm.vf.scale.PCA.cluster.umap.markers$avg_logFC)*(-log(s1f.sub.norm.vf.scale.PCA.cluster.umap.markers$p_val_adj)) ## rank for GSEA analysis later
write.csv(s1f.sub.norm.vf.scale.PCA.cluster.umap.markers,"s1f_clusters_DEG_genes.csv")


#this code works
FeaturePlot(s1f.sub.norm.vf.scale.PCA.cluster.umap, features = c("Zbtb20", "Map1b", "Erc2", "Olfm17", "Epha7", "Ppfia2", "Ppp3ca", "Rfx3",  "Adcy1"))
FeaturePlot(s1f.sub.norm.vf.scale.PCA.cluster.umap, features = c("Mbp"))
DotPlot(object = s1f.sub.norm.vf.scale.PCA.cluster.umap, features =  c("Mbp","Gria4","Ndst4","Rnf182"),cols = c("blue", "red")) + RotatedAxis()
VlnPlot(object = s1f.sub.norm.vf.scale.PCA.cluster.umap, features =  c("Mbp","Gria4","Ndst4","Rnf182"))

#sort by topl10 and dotplot it
top10 <- s1f.sub.norm.vf.scale.PCA.cluster.umap %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()
DotPlot(
  s1f.sub.norm.vf.scale.PCA.cluster.umap,
  assay = NULL,
  unique(top10$gene),
  cols = c("blue", "red"),
  col.min = 2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by= NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme (axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

#saveRDS(pbmc, file = "pbmc3k_final.rds")
#tools to help annotate
#http://bio-bigdata.hrbmu.edu.cn/CellMarker/search.jsp?cellMarkerSpeciesType-Human&cellMarker-pf4
#https://azimuth.hubmap.consortium.org/

cluster6 <- FindMarkers(s1f.sub.norm.vf.scale.PCA.cluster.umap, ident.1 = 6, ident.2 = c(4, 2,0), min.pct = 0.25)
cluster6$gene <-rownames (cluster6)

#visualize degs
EnhancedVolcano(cluster6,
                 lab = rownames(cluster6),
                 x = 'avg_log2FC',
                 y = 'p_val_adj',
                 title= 'Cluster 6 vs Cluster 4,2,0',
                 pCutoff = 10e-32,
                 FCcutoff = 1,
                 pointSize = 3.0,
                 labSize = 6.0)


# set up list of canonical cell type markers
canonical_markers <- list(
  'Astrocyte' = c('GFAP', 'AQP4', 'SLC1A2'),
  'Pan-neuronal' = c('SNAP25', 'SYT1'),
  'Excitatory Neuron' = c('SLC17A7', 'SATB2'),
  'Inhibitory Neuron' = c('GAD1', 'GAD2'),
  'Microglia' = c('CSF1R', 'CD74', 'P2RY12'),
  'Oligodendrocyte' = c('MOBP', 'MBP', 'MOG'),
  'Olig. Progenitor' = c('PDGFRA', 'CSPG4')
)

# plot heatmap:
library(viridis)
png('basic_canonical_marker_heatmap.png', width=10, height=10, units='in', res=200)
DoHeatmap(s1f.sub.norm.vf.scale.PCA.cluster.umap, group.by ="seurat_clusters", features=as.character(unlist(canonical_markers))) + scale_fill_gradientn(colors=viridis(256))
dev.off()
