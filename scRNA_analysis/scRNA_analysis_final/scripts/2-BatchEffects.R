
################################# START ########################################

## Load libraries --------------------------------------------------------------

library(pvca)
library(golubEsets)
library(Biobase)
library(Seurat)

## Load, Filter, Normalize Seurat object ---------------------------------------

seurat_object <- readRDS(obj.path)
counts <- as.matrix(GetAssayData(object = seurat_object, assay = "RNA", slot = "counts"))
seurat_object <- seurat_object[!rowSums(counts==0) > 0.7*ncol(counts),]
dim(seurat_object)
seurat_object <- NormalizeData(seurat_object)

## Create Expression set -------------------------------------------------------

# Seed
set.seed(1)
idx <- sample(seq(1,ncol(seurat_object),1),round(ncol(seurat_object)/10))
seurat_object_part <- seurat_object[,idx]
  
# Make expression set
exp <- as.matrix(GetAssayData(object = seurat_object_part,
                              assay = "RNA", slot = "data"))
pData=as.data.frame(seurat_object_part@meta.data)[,c("sample","donor","timepoint",
                                                     "treatment","Phase")]
metadata <- data.frame(labelDescription=c("Library Name","Patient ID","Treatment Timepoint","Ibrutinib Treatment",
                                          "Cell Cycle Phase"),row.names=c("sample","donor","timepoint","treatment","Phase"))
phenoData <- new("AnnotatedDataFrame",data=pData,
                 varMetadata=metadata)
ExpSet <- ExpressionSet(assayData=exp,
                            phenoData=phenoData)

## Create Expression set -------------------------------------------------------

pct_threshold <- 0.6
batch.list <- list()
batch.list[[1]] <- c("sample","treatment")
batch.list[[2]] <- c("donor","treatment")
batch.list[[3]] <- c("sample","treatment","Phase")
batch.list[[4]] <- c("timepoint", "donor","treatment","Phase")
batch.list[[5]] <- c("sample", "timepoint", "donor","treatment","Phase")

pdf(paste0(batch.dir,"batch_effects.pdf"))
for(i in 1:length(batch.list)){
  pvcaObj <- pvcaBatchAssess(ExpSet, batch.list[[i]], pct_threshold)
  bp <- barplot(pvcaObj$dat,  xlab = "Effects",
                ylab = "Weighted average proportion variance", ylim= c(0,1.1),
                col = c("blue"), las=2, main="PVCA estimation bar chart")
  axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
  values = pvcaObj$dat
  new_values = round(values , 3)
  text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)
}
dev.off()

################################## END #########################################

