#26.6.23 hippo samples analysis
library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
#load the >>>>>>>>dataset_1h
data1 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/1H_filtered_feature_bc_matrix")
H1 <- CreateSeuratObject(counts = data1,project = 'hippo', min.cells =3, min.features=200)

#load the >>>>>>>>dataset_2h
data2 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/2H_filtered_feature_bc_matrix")
H2 <- CreateSeuratObject(counts = data2,project = 'hippo', min.cells =3, min.features=200)

#load the >>>>>>>>dataset_3h
data3 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/3H_filtered_feature_bc_matrix")
H3 <- CreateSeuratObject(counts = data3,project = 'hippo', min.cells =3, min.features=200)


# QC
# % MT reads
H1[["percent.mt"]] <- PercentageFeatureSet(H1, pattern = "^Mt-")
hp1 <- VlnPlot(H1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(H1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(H1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


H2[["percent.mt"]] <- PercentageFeatureSet(H2, pattern = "^Mt-")
hp2 <- VlnPlot(H2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot3 <- FeatureScatter(H2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(H2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3, plot4))

H3[["percent.mt"]] <- PercentageFeatureSet(H3, pattern = "^Mt-")
hp3 <- VlnPlot(H3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hp3

plot5 <- FeatureScatter(H3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(H3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot5, plot6))

hp1+hp2+hp3

#Filtering
H1f <- subset(H1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
H2f <- subset(H2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
H3f <- subset(H3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

#merged the surat object
merged.H123 <- merge(H1f, y = c(H2f,H3f), add.cell.ids = c("1H", "2H", "3H"), project = "hippo123")

# head(colnames(merged.H123))
# tail(colnames(merged.H123))
# unique(sapply(X = strsplit(colnames(merged.H123), split = "_"), FUN = "[", 1))
# table(merged.H123$orig.ident)

#adding meta data
library(stringr)
sample <- names(merged.H123@active.ident)
sample_detect <- ifelse(str_detect(sample, "1H"), "1H", ifelse(str_detect(sample, "2H"), "2H", "3H"))

merged.H123@meta.data$sample <- sample_detect

Idents(object = merged.H123) <- "sample"

#Normalized Data and Intergrate the three dataset
mergedH123.list <-  SplitObject(merged.H123, split.by = "sample")

# perform standard preprocessing on each object
for (i in 1:length(mergedH123.list)) {
  mergedH123.list[[i]] <- NormalizeData(mergedH123.list[[i]], verbose = FALSE)
  mergedH123.list[[i]] <- FindVariableFeatures(
    mergedH123.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
  )
}

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = mergedH123.list)
mergedH123.list <- lapply(X = mergedH123.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Integrate Data
# find anchors
anchors <- FindIntegrationAnchors(object.list = mergedH123.list)
# integrate data
mergedH123.list <- IntegrateData(anchorset = anchors)

#Data Cleaning and Dimensional reduction for visualization
mergedH123.list <- ScaleData(mergedH123.list, verbose = FALSE)
mergedH123.list <- FindVariableFeatures(mergedH123.list, 
                                       selection.method = "vst",
                                       nfeatures = 2000, 
                                       verbose = FALSE)

mergedH123.list <- RunPCA(mergedH123.list, npcs = 30, verbose = FALSE)
mergedH123.list <- RunUMAP(mergedH123.list, reduction = "pca", dims = 1:30)
mergedH123.list <- FindNeighbors(mergedH123.list, reduction = "pca", dims = 1:30)
mergedH123.list <- FindClusters(mergedH123.list, resolution = 0.5)

# Visualization

p1 <- DimPlot(mergedH123.list, reduction = "umap",label = TRUE)
p1
Idents(object = mergedH123.list) <- "sample"
p2 <- DimPlot(mergedH123.list, reduction = "umap", cols = c("#cc0000", "#38761d","#9fc5e8"))
p2
p1+p2

#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
DimPlot(mergedH123.list, reduction = "umap", split.by = "sample")
