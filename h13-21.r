#27.6.23 zia
library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
#load the >>>>>>>>dataset_13h
data13 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/13H_filtered_feature_bc_matrix")
H13 <- CreateSeuratObject(counts = data13,project = 'hippo', min.cells =3, min.features=200)

#load the >>>>>>>>dataset_14h
data14 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/14H_filtered_feature_bc_matrix")
H14 <- CreateSeuratObject(counts = data14,project = 'hippo', min.cells =3, min.features=200)

#load the >>>>>>>>dataset_15h
data15 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/15H_filtered_feature_bc_matrix")
H15 <- CreateSeuratObject(counts = data15,project = 'hippo', min.cells =3, min.features=200)

#-----------------------------------------------------------------
#load the >>>>>>>>dataset_19h
data19 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/19H_filtered_feature_bc_matrix")
H19 <- CreateSeuratObject(counts = data15,project = 'hippo', min.cells =3, min.features=200)

#load the >>>>>>>>dataset_20h
data20 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/20H_filtered_feature_bc_matrix")
H20 <- CreateSeuratObject(counts = data15,project = 'hippo', min.cells =3, min.features=200)

#load the >>>>>>>>dataset_21h
data21 <- Read10X("C:/Users/ziaha/OneDrive/Documents/MobaXterm/home/21H_filtered_feature_bc_matrix")
H21 <- CreateSeuratObject(counts = data15,project = 'hippo', min.cells =3, min.features=200)


# QC
# % MT reads
H13[["percent.mt"]] <- PercentageFeatureSet(H13, pattern = "^Mt-")
hp13 <- VlnPlot(H13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hp13
plot1 <- FeatureScatter(H13, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(H13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


H14[["percent.mt"]] <- PercentageFeatureSet(H14, pattern = "^Mt-")
hp14 <- VlnPlot(H14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hp14
plot3 <- FeatureScatter(H14, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(H14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3, plot4))

H15[["percent.mt"]] <- PercentageFeatureSet(H15, pattern = "^Mt-")
hp15 <- VlnPlot(H15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hp15

plot5 <- FeatureScatter(H15, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(H15, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot5, plot6))
#------------------------------------------
H19[["percent.mt"]] <- PercentageFeatureSet(H19, pattern = "^Mt-")
hp19 <- VlnPlot(H19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hp19
plot7 <- FeatureScatter(H19, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot8 <- FeatureScatter(H19, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot7, plot8))

H20[["percent.mt"]] <- PercentageFeatureSet(H20, pattern = "^Mt-")
hp20 <- VlnPlot(H20, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hp20
plot9 <- FeatureScatter(H13, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot10 <- FeatureScatter(H13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot9, plot10))

H21[["percent.mt"]] <- PercentageFeatureSet(H21, pattern = "^Mt-")
hp21 <- VlnPlot(H21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hp21
plot11 <- FeatureScatter(H13, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot12 <- FeatureScatter(H13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot11, plot12))

hp1+hp2+hp3

#Filtering
H13f <- subset(H13, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
H14f <- subset(H14, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
H15f <- subset(H15, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)

#--------------------------------------------------------------------------------------------
H19f <- subset(H19, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
H20f <- subset(H20, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
H21f <- subset(H21, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)

#merged the surat object
merged.H131415 <- merge(H13f, y = c(H14f,H15f), add.cell.ids = c("13H", "14H", "15H"), project = "hippo131415")
merged.H192021 <- merge(H19f, y = c(H20f,H21f), add.cell.ids = c("19H", "20H", "21H"), project = "hippo192021")

merge.13_21 <- merge(merged.H131415, y = c(merged.H192021), add.cell.ids = c("cE14F", "cGFPF"), project = "hippo13_21")

# head(colnames(merge.13_21))
# tail(colnames(merge.13_21))
# unique(sapply(X = strsplit(colnames(merge.13_21), split = "_"), FUN = "[", 1))
# table($orig.ident)

#adding meta data
library(stringr)
sample <- names(merged.H131415@active.ident)
sample_detect <- ifelse(str_detect(sample, "13H"), "13H", ifelse(str_detect(sample, "14H"), "14H", "15H"))

sample <- names(merged.H192021@active.ident)
sample_detect <- ifelse(str_detect(sample, "19H"), "19H", ifelse(str_detect(sample, "20H"), "20H", "21H"))

# for two mega merged group
sample <- names(merge.13_21@active.ident)
sample_detect <- ifelse(str_detect(sample, "cE14F"), "cE14F", "cGFPF")

merge.13_21@meta.data$sample <- sample_detect
Idents(object = merge.13_21) <- "sample"

merged.H131415@meta.data$sample <- sample_detect
Idents(object = merged.H131415) <- "sample"

merged.H192021@meta.data$sample <- sample_detect
Idents(object = merged.H192021) <- "sample"

#Normalized Data and Intergrate the three dataset
merged.H131415.list <-  SplitObject(merged.H131415, split.by = "sample")
merged.H192021.list <-  SplitObject(merged.H192021, split.by = "sample")
merge.13_21.list <-  SplitObject(merge.13_21, split.by = "sample")


# perform standard preprocessing on each object
for (i in 1:length(merge.13_21.list)) {
  merge.13_21.list[[i]] <- NormalizeData(merge.13_21.list[[i]], verbose = FALSE)
  merge.13_21.list[[i]] <- FindVariableFeatures(
    merge.13_21.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
  )
}

for (i in 1:length(merged.list)) {
  merged.H131415.list[[i]] <- NormalizeData(merged.H131415.list[[i]], verbose = FALSE)
  merged.H131415.list[[i]] <- FindVariableFeatures(
    merged.H131415.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
  )
}

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = merged.H131415.list)
merged.H131415.list <- lapply(X = merged.H131415.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


features <- SelectIntegrationFeatures(object.list = merge.13_21.list)
merge.13_21.list <- lapply(X = merge.13_21.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Integrate Data
# find anchors
anchors <- FindIntegrationAnchors(object.list = merged.H131415.list)
anchors <- FindIntegrationAnchors(object.list = merge.13_21.list)

# integrate data
merged.H131415.list <- IntegrateData(anchorset = anchors)

merge.13_21.list <- IntegrateData(anchorset = anchors)

#Data Cleaning and Dimensional reduction for visualization
merged.H131415.list <- ScaleData(merged.H131415.list, verbose = FALSE)
merged.H131415.list <- FindVariableFeatures(merged.H131415.list, 
                                       selection.method = "vst",
                                       nfeatures = 2000, 
                                       verbose = FALSE)

merged.H131415.list <- RunPCA(merged.H131415.list, npcs = 30, verbose = FALSE)
merged.H131415.list <- RunUMAP(merged.H131415.list, reduction = "pca", dims = 1:30)
merged.H131415.list <- FindNeighbors(merged.H131415.list, reduction = "pca", dims = 1:30)
merged.H131415.list <- FindClusters(merged.H131415.list, resolution = 0.5)
gc()

# for two mega merged group--------------------------------------------
#Data Cleaning and Dimensional reduction for visualization
merged.H131415.list <- ScaleData(merged.H131415.list, verbose = FALSE)
merged.H131415.list <- FindVariableFeatures(merged.H131415.list, 
                                            selection.method = "vst",
                                            nfeatures = 2000, 
                                            verbose = FALSE)

merged.H131415.list <- RunPCA(merged.H131415.list, npcs = 30, verbose = FALSE)
merged.H131415.list <- RunUMAP(merged.H131415.list, reduction = "pca", dims = 1:30)
merged.H131415.list <- FindNeighbors(merged.H131415.list, reduction = "pca", dims = 1:30)
merged.H131415.list <- FindClusters(merged.H131415.list, resolution = 0.5)
gc()



# Visualization

p1 <- DimPlot(merged.H131415.list, reduction = "umap",label = TRUE)
p1
Idents(object = merged.H131415.list) <- "sample"
p2 <- DimPlot(merged.H131415.list, reduction = "umap", cols = c("#cc0000", "#38761d","#9fc5e8"))
p2
p1+p2

#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
DimPlot(merged.H131415.list, reduction = "umap", split.by = "sample")
