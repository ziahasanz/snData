# Zia - April 13th, 2022

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


s1f.sub <- subset(s1f, subset = nFeature_RNA > 100 & nCount_RNA > 600 & percent.Mt < 5)
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


#Determine the ‘dimensionality’ of the dataset
#if you use below two commands for the‘dimensionality’, it will take longtime to analyze the data (for instance: 4h for my data). I skip this step and use only the Elbow plot.
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
cluster0.markers <- FindMarkers(object = s1f.sub.norm.vf.scale.PCA.cluster.tSNE, ident.1 = 1, min.pct = 0.25)
head(x = cluster0.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
s1f.sub.norm.vf.scale.PCA.cluster.tSNE.markers <- FindAllMarkers(object = s1f.sub.norm.vf.scale.PCA.cluster.tSNE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(s1f.sub.norm.vf.scale.PCA.cluster.tSNE.markers)
s1f.sub.norm.vf.scale.PCA.cluster.tSNE.markers %>%group_by(cluster) %>%slice_max(n = 2, order_by = avg_log2FC)
VlnPlot(object = s1f.sub.norm.vf.scale.PCA.cluster.tSNE, features = c("RBFOX3"))
