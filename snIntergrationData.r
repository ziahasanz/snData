#zia hasan 20th July
setwd("single_integration")
library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)

#load the >>>>>>>>dataset1
data1 <- Read10X("../single_integration/raw_feature_bc_matrix_s1f_cE14")
s1f.ce14 <- CreateSeuratObject(counts = data1,project = 'hippo', min.cells =3, min.features=200)
# QC
# % MT reads
s1f.ce14[["percent.mt"]] <- PercentageFeatureSet(s1f.ce14, pattern = "^Mt-")
VlnPlot(s1f.ce14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(s1f.ce14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
#Filtering
s1f.ce14 <- subset(s1f.ce14, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#load the >>>>>>>dataset2
data2 <- Read10X("../single_integration/raw_feature_bc_matrix_s2f_cGFP")
s2f.cgfp <- CreateSeuratObject(counts = data2,project = 'hippo', min.cells =3, min.features=200)
# QC
# % MT reads
s2f.cgfp[["percent.mt"]] <- PercentageFeatureSet(s2f.cgfp, pattern = "^Mt-")
VlnPlot(s2f.cgfp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Filtering
s2f.cgfp <- subset(s2f.cgfp, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#merged the surat object
merged.ce14gfp <- merge(x = s1f.ce14, y = s2f.cgfp, add.cell.ids = c("cE14", "cGFP"))  

#adding meta data
library(stringr)
sample <- names(merged.ce14gfp@active.ident)
sample_detect <- ifelse(str_detect(sample,"cE14"),"cE14","cGFP")

merged.ce14gfp@meta.data$sample <- sample_detect

Idents(object = merged.ce14gfp) <- "sample"

#Normalized Data and Intergrate the two dataset
ce14gfp.list <-  SplitObject(merged.ce14gfp, split.by = "sample")
# perform standard preprocessing on each object
for (i in 1:length(ce14gfp.list)) {
  ce14gfp.list[[i]] <- NormalizeData(ce14gfp.list[[i]], verbose = FALSE)
  ce14gfp.list[[i]] <- subset(ce14gfp.list[[i]], downsample = 10000)
  ce14gfp.list[[i]] <- FindVariableFeatures(
    ce14gfp.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
  )
}


# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ce14gfp.list)
ce14gfp.list <- lapply(X = ce14gfp.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


#Intergrate Data
# find anchors
anchors <- FindIntegrationAnchors(object.list = ce14gfp.list)
# integrate data
merged.ce14gfp <- IntegrateData(anchorset = anchors)

#Data Cleaning and Dimensional reduction for visualization
merged.ce14gfp <- ScaleData(merged.ce14gfp, verbose = FALSE)
merged.ce14gfp <- FindVariableFeatures(merged.ce14gfp, 
                                       selection.method = "vst",
                                       nfeatures = 2000, 
                                       verbose = FALSE)

merged.ce14gfp <- RunPCA(merged.ce14gfp, npcs = 30, verbose = FALSE)
merged.ce14gfp <- RunUMAP(merged.ce14gfp, reduction = "pca", dims = 1:30)
merged.ce14gfp <- FindNeighbors(merged.ce14gfp, reduction = "pca", dims = 1:30)
merged.ce14gfp <- FindClusters(merged.ce14gfp, resolution = 0.5)
gc()

# Visualization

p1 <- DimPlot(merged.ce14gfp, reduction = "umap",label = TRUE)
p1

Idents(object = merged.ce14gfp) <- "sample"
p2 <- DimPlot(merged.ce14gfp, reduction = "umap", cols = c("red", "green"))
p2
p1+p2

#nucelei count in each clusters
merged.ce14gfp.cellnos <- Seurat::StashIdent(object = merged.ce14gfp, save.name = "clusters.cellnos")
table(merged.ce14gfp.cellnos@meta.data$clusters.cellnos, merged.ce14gfp.cellnos@meta.data$orig.ident)

# Find differentially expressed features between ....
# search for positive markers
Idents(object = merged.ce14gfp) <- "sample"
sammple.markers <- FindMarkers(merged.ce14gfp, ident.1 = "cE14", ident.2 = "cGFP")
# view results
head(sammple.markers)

# DEGs markers
library(ComplexHeatmap)
library(BiocManager)

heatmapdf <- sammple.markers[1:25,]
row_ha = rowAnnotation("CE14" = anno_barplot(heatmapdf$pct.1),
                       "cGFP"= anno_barplot(heatmapdf$pct.2),
                       width = unit(10, "cm"))

ht0 <- Heatmap(heatmapdf$avg_log2FC,
               name = "Log2FC",
               cluster_rows = TRUE, 
               row_labels = rownames(heatmapdf), 
               right_annotation = row_ha,
               width = unit(1, "cm"))

ht0

#Find differentially expressed features between one group and  other group
#search for positive markers
Idents(object = merged.ce14gfp) <- "seurat_clusters"
cluster0.de.markers <- FindMarkers(merged.ce14gfp, ident.1 = "0", ident.2 = NULL, only.pos = TRUE)
head(cluster0.de.markers)

cluster1.de.markers <- FindMarkers(merged.ce14gfp, ident.1 = "1", ident.2 = NULL, only.pos = TRUE)
head(cluster1.de.markers)

cluster01.de.markers <- FindMarkers(merged.ce14gfp, ident.1 = "0", ident.2 = "1", only.pos = TRUE)
head(cluster01.de.markers)

#Identify conserved cell type markers (when you have two different samples. EX: control and Experimental group)
DefaultAssay(merged.ce14gfp) <- "RNA"

cluster0.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 0, grouping.var = "sample", verbose = FALSE)
head(cluster0.markers)

cluster1.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 1, grouping.var = "sample", verbose = FALSE)
head(cluster1.markers)

cluster2.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 2, grouping.var = "sample", verbose = FALSE)
head(cluster2.markers)


cluster3.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 3, grouping.var = "sample", verbose = FALSE)
head(cluster3.markers)

cluster4.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 4, grouping.var = "sample", verbose = FALSE)
head(cluster4.markers)


cluster5.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 5, grouping.var = "sample", verbose = FALSE)
head(cluster5.markers)

cluster6.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 6, grouping.var = "sample", verbose = FALSE)
head(cluster6.markers)

cluster7.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 7, grouping.var = "sample", verbose = FALSE)
head(cluster7.markers)

cluster8.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 8, grouping.var = "sample", verbose = FALSE)
head(cluster8.markers)

cluster9.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 9, grouping.var = "sample", verbose = FALSE)
head(cluster9.markers)

cluster10.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 10, grouping.var = "sample", verbose = FALSE)
head(cluster10.markers)

cluster11.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 11, grouping.var = "sample", verbose = FALSE)
head(cluster11.markers)

cluster12.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 12, grouping.var = "sample", verbose = FALSE)
head(cluster12.markers)

cluster13.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 13, grouping.var = "sample", verbose = FALSE)
head(cluster13.markers)

cluster14.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 14, grouping.var = "sample", verbose = FALSE)
head(cluster14.markers)

cluster15.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 15, grouping.var = "sample", verbose = FALSE)
head(cluster15.markers)

cluster16.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 16, grouping.var = "sample", verbose = FALSE)
head(cluster16.markers)

cluster17.markers <- FindConservedMarkers(merged.ce14gfp, ident.1 = 17, grouping.var = "sample", verbose = FALSE)
head(cluster17.markers)


#Cluster identifiaction using known markers

##Cholinergic neurons
VlnPlot(object = merged.ce14gfp, features =  c("Chat","Slc18a3","Ache"), pt.size = 0)
#Dopaminergic neurons
VlnPlot(object = merged.ce14gfp, features =  c("Th","Slc6a3","Foxa2","Kcnj6","Nr4a2","Lmx1b","Dbh","Slc6a2","Ppp1r1b"), pt.size = 0)
#Serotonergic neurons
VlnPlot(object = merged.ce14gfp, features =  c("Tph1","Slc6a4","Fev","Htr1a","Htr1b"), pt.size = 0)
#GABAergic neurons 8,9
VlnPlot(object = merged.ce14gfp, features =  c("Slc6a1","Gabbr1","Gabbr2","Gad2","Gad1","Slc32a1"), pt.size = 0)
#Glutamatergic neurons 3,5,7,10,11,13,14,15?
VlnPlot(object = merged.ce14gfp, features =  c("Slc17a7","Slc17a6","Grin1","Grin2b","Gls","Glul","Grin2a"), pt.size = 0)
#Only for Neurons
VlnPlot(object = merged.ce14gfp, features =  c("Camk2a","Kif5c"), pt.size = 0)

#Mature neurons 15?
VlnPlot(object = merged.ce14gfp, features =  c("Rbfox3","Map2","Syp","Dlg4","Tubb3","Mapt","Ina","Gap43","Nrp1"), pt.size = 0)
#Microglia 6
VlnPlot(object = merged.ce14gfp, features =  c("P2ry12","Itgam","Cd40","Ptprc","Cd68","Aif1","Cx3cr1","Tmem119","Adgre1","C1qa","Nos2","Tnf","Isyna1","Ccl4","Adora3","Adrb2","Bhlhe41","Bin1","Klf2","Nav3","Rhob","Sall1","Siglec8","Slc1a3","Spry1","Tal1"), pt.size = 0)
#Astrocytes 12
VlnPlot(object = merged.ce14gfp, features =  c("Gfap", "Slc1a3", "Slc1a2", "Glul", "S100b", "Aldh1l1", "Aqp4", "Igfbp3", "Atp13a4", "Cbs", "Sox9", "Cd40", "Cd80", "Mfge8", "Cd86", "C5ar1"), pt.size = 0)
#Oligodendrocytes 17,2,
VlnPlot(object = merged.ce14gfp, features =  c("Ng2","Olig1","Olig2","Olig3","Cldn11","Mbp","Mog","Mag","Galc","Cnp","Sox10","Fa2h","Ugt8"), pt.size = 0)
#Oligodendrocytes precursor cells 4
VlnPlot(object = merged.ce14gfp, features =  c("Lhfpl3","Megf11","Pcdh15","Pdgfra","Cspg4","Rnf5"), pt.size = 0)

---------------------------------------------------------------------
#Further confirmation of the cluster identity using top markers (FindConservedmarker) from each cluster (based on logfold change, pct.1 and pct2)

#Cluster 0 Glutamatergic neurons?
VlnPlot(object = merged.ce14gfp, features = c("Ptgds","Hsp90ab1","Atp1b1","Actb","Hsp90ab1"), pt.size = 0.5)

DotPlot(object = merged.ce14gfp, features = c("Ptgds","Hsp90ab1","Atp1b1","Actb"), cols = c("blue", "red")) + RotatedAxis()

#Cluster 1 GABAergic neurons?
VlnPlot(object = merged.ce14gfp, features =  c("Ubb","Tmsb10","Ptgds","Plp1","Apoe","Atp6v0a1","Calm1.1"), pt.size = 0.5)
FeaturePlot(object = merged.ce14gfp, features = c("Ptgds","Ldhb","Rnasek","Ptms"), min.cutoff = "q10", pt.size = 0.5, label = TRUE, cols = c("grey", "red","darkred"))

#Cluster 2 Oligodendrocytes
VlnPlot(object = merged.ce14gfp, features =  c("Luzp2","Sorcs1","Pde8a","Map7","Nkain2","Qk","Fnbp1","Mast4","Grm3"), pt.size = 0)
FeaturePlot(object = merged.ce14gfp, features = c("Luzp2","Sorcs1","Pde8a","Map7","Nkain2","Qk","Fnbp1","Grm3","Mast4"), min.cutoff = "q10", pt.size = 0.5, cols = c("grey", "red","darkred"))

#cluster 3 Glutamatergic neurons (inhibitory)
VlnPlot(object = merged.ce14gfp, features =  c("Luzp2","Sorcs1","Pde8a","Map7","Nkain2","Qk","Fnbp1","Mast4","Grm3"), pt.size = 0)
#cluster 4 Oligodendrocytes precursor cells
#cluster 5 Glutamatergic neurons (inhibitory?)
#cluster 6 Microglia
#cluster 7 Glutamatergic Neuron (Excitatory)
#cluster 8 GABAergic neurons (inhibitory)
#cluster 9 GABAergic neurons (inhibitory)
#cluster 10 Glutamatergic Neurons (Excitatory)
#cluster 11 Glutamatergic Neurons(Excitatory)
VlnPlot(object = merged.ce14gfp, features =  c("Tshz2","Dpp10","Hs6st3","Sgcz","Pdzrn3","Vwc2l","Satb2","Tox","Fhod3"), pt.size = 0)
FeaturePlot(object = merged.ce14gfp, features = c("Tshz2","Dpp10","Hs6st3","Sgcz","Pdzrn3","Vwc2","Satb2","Tox","Fhod3"), min.cutoff = "q10", pt.size = 0.5, cols = c("grey", "red","darkred"))
DotPlot(object = merged.ce14gfp, features =  c("Tshz2","Dpp10","Hs6st3","Sgcz","Pdzrn3","Vwc2","Satb2","Tox","Fhod3"),cols = c("blue", "red")) + RotatedAxis()
#Gaba - Hs6st3, Gluta-Vwc2l,Tshz2

#cluster 12 Astrocytes
#cluster 13 Glutamatergic Neurons(Excitatory)
#cluster 14 GABAergic neurons
#cluster 15 GABAergic neurons
#cluster 16 Pericytes
VlnPlot(object = merged.ce14gfp, features =  c("Utrn", "Ssh2","Dlc1","Abcc9","Atp13a5","Ebf1","Egflam","Rgs5","Mecom"), pt.size = 0)

#cluster 17 Oligodendrocytes

-------------------------------------------------
#Rename the clusters
merged.ce14gfp.renamed <- RenameIdents(merged.ce14gfp, "0" = "Glutamatergic1", "1" = "GABAergic1", "2" = "Oligodendrocytes1","3" = "Glutamatergic2", "4" = "OPC", "5" = "Glutamatergic3", "6" = "Microglia", "7" = "Glutamatergic4", "8" = "GABAergic2", "9" = "GABAergic3", "10" = "Glutamatergic5", "11" = "Glutamatergic6", "12" = "Astrocytes", "13" = "Glutamatergic7", "14" = "GABAergic4", "15" = "GABAergic5", "16" = "Pericytes", "17" = "Oligodendrocytes2")
DimPlot(merged.ce14gfp.renamed, label = TRUE, label.size = 3)

#viewing conserved cell type markers across conditions
Idents(merged.ce14gfp.renamed) <- factor(Idents(merged.ce14gfp.renamed), levels = c("Glutamatergic1", "GABAergic1", "Oligodendrocytes1", "Glutamatergic2", "Oligodendrocytes2",  "Glutamatergic3",  "Microglia",  "Glutamatergic4",  "GABAergic2",  "GABAergic3","Glutamatergic5",  "Glutamatergic6",  "Astrocytes",  "Glutamatergic7",  "GABAergic4", "GABAergic5",   "Pericytes",  "Oligodendrocytes3"))

#GABA
markers.to.plot.gaba <- c("Slc6a1","Gabbr1","Gabbr2","Gad2","Gad1","Slc32a1")
DotPlot(merged.ce14gfp.renamed, features = markers.to.plot.gaba, cols = c("blue", "red"), split.by = "sample") + RotatedAxis()


#Glutamate
markers.to.plot.glu <- c("Zbtb20","Tshz2","Dpp10","Hs6st3","Sgcz","Pdzrn3","Vwc2l","Satb2","Tox","Fhod3")
DotPlot(merged.ce14gfp.renamed, features = markers.to.plot.glu, cols = c("blue", "red"), split.by = "sample") + RotatedAxis()

#Astrocytes
markers.to.plot.astro <- c("Gfap", "Slc1a3", "Slc1a2", "Glul", "S100b", "Aldh1l1", "Aqp4", "Igfbp3", "Atp13a4", "Cbs", "Sox9", "Cd40", "Cd80", "Mfge8", "Cd86", "C5ar1")
DotPlot(merged.ce14gfp.renamed, features = markers.to.plot.astro, cols = c("blue", "red"), split.by = "sample") + RotatedAxis()

#Microglia
markers.to.plot.micro <- c("P2ry12","Ptprc","Srgap2","Cx3cr1","Arhgap15","Tmcc3","Inpp5d","Zfhx3","Slc1a3","Blnk","Slco2b1","Dock8","Cped1","Mertk")
DotPlot(merged.ce14gfp.renamed, features = markers.to.plot.micro, cols = c("blue", "red"), split.by = "sample") + RotatedAxis()

#Oligodendrocytes
markers.to.plot.oligo <- c("Lhfpl3","Megf11","Pcdh15","Pdgfra","Cspg4","Rnf5")
DotPlot(merged.ce14gfp.renamed, features = markers.to.plot.oligo, cols = c("blue", "red"), split.by = "sample") + RotatedAxis()


#OPC
markers.to.plot.oligo <- c("Lhfpl3","Megf11","Pcdh15","Pdgfra","Cspg4","Rnf5")
DotPlot(merged.ce14gfp.renamed, features = markers.to.plot.oligo, cols = c("blue", "red"), split.by = "sample") + RotatedAxis()
VlnPlot(merged.ce14gfp.renamed, features = markers.to.plot.oligo, cols = c("blue", "red"), split.by = "sample", pt.size = 0)

#Pericytes
markers.to.plot.peri<- c("Utrn", "Ssh2","Dlc1","Abcc9","Atp13a5","Ebf1","Egflam","Rgs5","Mecom")
DotPlot(merged.ce14gfp.renamed, features = markers.to.plot.peri, cols = c("blue", "red"), split.by = "sample") + RotatedAxis()

#Identify differential expressed genes across conditions
