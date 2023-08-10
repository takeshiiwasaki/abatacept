library(dplyr)
library(Seurat)
library(patchwork)
pbmc <- as.matrix(read.table(file = "~/AMP/RA/exprMatrix.tsv" ,row.names = 1 , header = TRUE))
pbmc <- CreateSeuratObject(counts = pbmc , project = "RA")

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
## Computing nearest neighbor graph
## Computing SNN
pbmc <- FindClusters(pbmc, resolution = 0.5)

#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap",label=TRUE)

#T
FeaturePlot(pbmc, features = c("CD3D"),label=TRUE)
#0,4,8:T細胞
#CD8+
FeaturePlot(pbmc, features = c("CD8A"),label=TRUE)
#4:CD8
#Naive C4+T
FeaturePlot(pbmc, features = c("IL7R"),label=TRUE)
FeaturePlot(pbmc, features = c("CCR7"),label=TRUE)
#8:Naive CD4+
#Memory CD4+
FeaturePlot(pbmc, features = c("IL7R"),label=TRUE)
FeaturePlot(pbmc, features = c("S100A4"),label=TRUE)
#0:Memory CD4+
#NK
FeaturePlot(pbmc, features = c("GNLY"),label=TRUE)
FeaturePlot(pbmc, features = c("NKG7"),label=TRUE)
#NK cell is inculded in cluster 4

#CD14+Mono/Macrophage, 
FeaturePlot(pbmc, features = c("CD14"),label=TRUE)
FeaturePlot(pbmc, features = c("LYZ"),label=TRUE)
#3:Classical monocyte 
#FCGR3A+Mono
FeaturePlot(pbmc, features = c("FCGR3A"),label=TRUE)
FeaturePlot(pbmc, features = c("MS4A7"),label=TRUE)
#9:Non classical monocte

#B
FeaturePlot(pbmc,features=c("MS4A1"),label=TRUE)
#2,10:B細胞

#DC
FeaturePlot(pbmc, features = c("FCER1A"),label=TRUE)
FeaturePlot(pbmc, features = c("CST3"),label=TRUE)

#Plasmablast marker (Nature Medicine 26,1070-1076 (2020))
FeaturePlot(pbmc,features=c("CD27"),label=TRUE)
FeaturePlot(pbmc,features=c("CD38"),label=TRUE)
FeaturePlot(pbmc,features=c("TNFRSF17"),label=TRUE)
#7:Plasmablast

#Fibroblast
FeaturePlot(pbmc,features=c("COL1A1"),label=TRUE)
FeaturePlot(pbmc,features=c("COL1A2"),label=TRUE)
FeaturePlot(pbmc,features=c("DCN"),label=TRUE)
#1,5,6:Fibroblast

#Platelet
#VlnPlot(pbmc, features = c("PPBP"))
#PPBP was not quantified in this dataset.

saveRDS(pbmc, file = "~/AMP/RA/result.rds")
pbmc <- RenameIdents(pbmc, `0` = "Mem CD4 T", `1` = "Fibroblast", `2` = "B",`3` = "CL Mono", `4` = "CD8 T", `5` = "Fibroblast", `6` = "Fibroblast",`7`="Plasmablast",`8`="Naive CD4 T",`9`="NC Mono",`10`="B")
