---
title: "Untitled"
output: html_document
---

```{r}
# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE68215", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL4133", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
par(mar=c(7,4,2,1))
title <- paste ("GSE68215", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE68215", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE68215")

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
write.table(ex,file="/mnt/c/Users/tiwasaki/Documents/RA_RNA.ABT/GSE68215/expr.xls",quote=FALSE,sep="\t")
pheno=as.data.frame(gset@phenoData@data)
write.table(pheno,file="/mnt/c/Users/tiwasaki/Documents/RA_RNA.ABT/GSE68215/pheno.xls",quote=FALSE,sep="\t")
```
```{r}
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "agilent_wholegenome_4x44k_v2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "agilent_wholegenome_4x44k_v2",
  values = rownames(exprs(gset)),
  uniqueRows=TRUE)
write.table(as.data.frame(annotLookup),file="/mnt/c/Users/tiwasaki/Documents/RA_RNA.ABT/GSE68215/annot.GSE68215.xls",quote=FALSE,sep="\t", col.names = NA,)
```
```{r}
#https://wiki.bits.vib.be/index.php/Analyse_GEO2R_data_with_R_and_Bioconductor
#G1:R, G0:NR
sml <- c("G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0");
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
 
# the differential expression contrast is set here to "NR" versus "R"
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))

write.table(tT, file="/mnt/c/Users/tiwasaki/Documents/RA_RNA.ABT/GSE68215/GSE68215.xls", 
            quote = FALSE, sep="\t", col.names = NA, row.names = T)

plot(tT$logFC, 1-tT$adj.P.Val, xlim=c(-6, 6), 
     main="Non-responder vs. responder",
     xlab="log2Ratio", ylab="1-adj.P.Val")
abline(h=0.95, col="red")
```
