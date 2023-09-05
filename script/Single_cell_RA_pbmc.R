library(Seurat)
library(SeuratData)
library(patchwork)
pbmc <- Read10X_h5("/home/tiwasaki/RA_pbmc/GSM4819747_RA_filtered_feature_bc_matrix.h5")

pbmc <- CreateSeuratObject(counts = pbmc, project = "RA", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##Cell QC
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
#Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#Dimention reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:20)

#T cell
VlnPlot(pbmc, features = c("CD3D"))
#0,1,2,3,9: T cells

#CD8+
VlnPlot(pbmc, features = c("CD8A"))
#0,1:CD8

#Naive C4+T
VlnPlot(pbmc, features = c("IL7R","CCR7"))
#Memory CD4+
VlnPlot(pbmc, features = c("IL7R","S100A4"))
#2,3,9:CD4+

#B
VlnPlot(pbmc, features = c("MS4A1"))
#6,7,12:B

#NK
VlnPlot(pbmc, features = c("GNLY","NKG7"))
#5:NK

#Mono
VlnPlot(pbmc, features = c("CD14","LYZ"))
#4,8:Mono

#Filter cluster with low number of cells
pbmc <- subset(pbmc, subset=seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9))

pbmc <- RenameIdents(pbmc, `0` = "CD8", `1` = "CD8", `2` = "CD4",`3` = "CD4", `4` = "Mono", `5` = "NK", `6` = "B", `7` = "B", `8` = "Mono", `9` = "CD4")
DimPlot(pbmc, label = TRUE)
table(Idents(pbmc))
saveRDS(pbmc, file = "RA_pbmc/GSM4819747_RA_filtered_feature_bc_matrix.rds")
