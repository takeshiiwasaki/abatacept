library(MetaVolcanoR) 
GSE78068=read.table("/mnt/c/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/limma_output.xls",header=1)
GSE68215=read.table("/mnt/c/Users/tiwasaki/Documents/RA_RNA.ABT/GSE68215/limma_output.xls",header=1)
GSE172188=read.table("/mnt/c/Users/tiwasaki/Documents/RA_RNA.ABT/GSE172188/limma_output.xls",header=1)
d=list()
d$GSE78068=GSE78068
d$GSE68215=GSE68215
d$GSE172188=GSE172188
#d$KURAMA=KURAMA
meta_degs_comb <- combining_mv(diffexp=d,
			       pcriteria='P.Value', 
			       foldchangecol='LogFC',
			       genenamecol='ID',
			       geneidcol=NULL,
			       metafc='Mean',
			       metathr=0.01, 
			       collaps=TRUE,
			       jobname="MetaVolcano",
			       outputfolder=".",
			       draw='HTML')

# Combining results
head(meta_degs_comb@metaresult, 3)
write.table(as.data.frame(meta_degs_comb@metaresult), file="/mnt/c/Users/tiwasaki/Documents/RA_RNA.ABT/meta/meta.xls", 
            quote = FALSE, sep="\t", col.names = NA, row.names = T)
