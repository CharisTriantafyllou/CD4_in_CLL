

################################# START ########################################

## Load libraries --------------------------------------------------------------

library(Seurat) # Seurat related utilities
library(SeuratData) # Seurat related utilities
library(SeuratObject) # Seurat related utilities
library(sctransform) # SCTransform code and utilities
library(glmGamPoi) # SCTransform model related code
library(SingleCellExperiment) # SCE object related utilities
library(SingleR) # Annotation package
library(celldex) # Annotation related utilities
library(ggplot2) # Plot related utilities
library(dplyr) # Data wrangling
library(RColorBrewer) # For advanced color Palette options
library(patchwork) # For combining plots
library(scater) # Optional single cell object utilities
library(scuttle) # Optional single cell object utilities

## Load the Seurat Object ------------------------------------------------------

# Messages for 'norm' mode
if(processing.mode=="norm"){
  message("Processing mode: normalization before annotation")
  prefix="norm"
  message("Loading postQC object . . .")
  seurat_path <- obj.path
  seurat_object <- readRDS(seurat_path)
  seurat_object@project.name <- paste0(prefix,"_Seurat")
  
  # Define outdir
  outdir=paste0(tr_cl_ann.dir,prefix,"/")
  dir.create(outdir,recursive = T)
  
}

# Messages for 'sct' mode
if(processing.mode=="sct"){
  message("Processing mode: sct before annotation")
  prefix="sct"
  message("Loading postQC object . . .")
  seurat_path <- obj.path
  seurat_object <- readRDS(seurat_path)
  seurat_object@project.name <- paste0(prefix,"_Seurat")
  
  #Transform
  message("Normalizing and transforming object . . .")
  seurat_object <- SCTransform(seurat_object, vst.flavor = "v2", verbose = FALSE)
  
  # Define outdir
  outdir=paste0(tr_cl_ann.dir,prefix,"/")
  dir.create(outdir,recursive = T)
  
}

# Messages for 'sct_sample_correction' mode
if(processing.mode=="sct_sample_correction"){
  message("Processing mode: normalization and sample correction before sct transformation and before annotation")
  prefix="sct_sample_correction"
  message("Loading postQC object . . .")
  seurat_path <- obj.path
  seurat_object <- readRDS(seurat_path)
  seurat_object@project.name <- paste0(prefix,"_Seurat")
  
  #Transform
  message("Normalizing and transforming object . . .")
  seurat_object <- SCTransform(seurat_object, vars.to.regress = c("sample"), vst.flavor = "v2", verbose = FALSE)
  
  # Define outdir
  outdir=paste0(tr_cl_ann.dir,prefix,"/")
  dir.create(outdir,recursive = T)
  
}

# Messages for 'sct_all_corrections' mode
if(processing.mode=="sct_all_corrections"){
  message("Processing mode: normalization and correction for all batch variables before sct transformation and before annotation")
  prefix="sct_all_corrections"
  message("Loading postQC object . . .")
  seurat_path <- obj.path
  seurat_object <- readRDS(seurat_path)
  seurat_object@project.name <- paste0(prefix,"_Seurat")
  
  #Transform
  message("Normalizing and transforming object . . .")
  seurat_object <- SCTransform(seurat_object, vars.to.regress = c("sample","donor","timepoint"), vst.flavor = "v2", verbose = FALSE)
  
  # Define outdir
  outdir=paste0(tr_cl_ann.dir,prefix,"/")
  dir.create(outdir,recursive = T)
  
}

# Messages for 'sct_integration' mode
if(processing.mode=="sct_integration"){
  message("Processing mode: integration after sct before annotation")
  prefix="sct_integration"
  message("Loading postQC object . . .")
  seurat_path <- obj.path
  seurat_object <- readRDS(seurat_path)
  seurat_object@project.name <- paste0(prefix,"_Seurat")
  
  # Define outdir
  outdir=paste0(tr_cl_ann.dir,prefix,"/")
  dir.create(outdir,recursive = T)
  
  # split the dataset into a list of two seurat objects (stim and CTRL)
  #object.list <- SplitObject(seurat_object, split.by = "sample")
  
  # Normalize object
  #message("Normalizing and transforming object . . .")
  
  # normalize and identify variable features for each dataset independently
  #object.list <- lapply(X = object.list, FUN = function(x) {
  #  x <- SCTransform(x, vst.flavor="v2",ncells = ncol(x), verbose = FALSE)
  #})
  
  # select features that are repeatedly variable across datasets for integration
  #features <- SelectIntegrationFeatures(object.list = object.list)
  #object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = features)
  #immune.anchors <- FindIntegrationAnchors(object.list = object.list, normalization.method = "SCT",
  #                                        anchor.features = features)
  #immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
  #saveRDS(immune.combined.sct,paste0(outdir,prefix,"_processed_seurat.rds"))
  immune.combined.sct <- readRDS(paste0(outdir,prefix,"_processed_seurat.rds"))
  seurat_object <- immune.combined.sct
  
  if(use_integrated_assay==T){
    message("Using integrated assay . . .")
    DefaultAssay(seurat_object) <- "integrated"
    prefix=paste0("integrated_",prefix)
  }else if(use_integrated_assay==F){
    message("Using SCT assay . . .")
    DefaultAssay(seurat_object) <- "SCT"
    prefix=paste0("SCT_",prefix)
  }
  
}


## Make plots ------------------------------------------------------------------


# Define references
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
bpe.ref <- celldex::BlueprintEncodeData()

# Create single cell experiments for annotation
sce <- as.SingleCellExperiment(seurat_object)

if(processing.mode=="norm"){
  sce_norm <- scuttle::logNormCounts(sce)
}else{
  sce_norm <- sce #make sure the "SCT assay is used though"
}

# Cell based annotation

hpca <- SingleR(sce_norm, ref=hpca.ref, labels=hpca.ref$label.main,
                 assay.type.test = "logcounts", assay.type.ref = "logcounts")
saveRDS(hpca, paste0(outdir,prefix,"_hpca_main_annotated.RDS"))

bpe <- SingleR(sce_norm, ref=bpe.ref, labels=bpe.ref$label.main,
                 assay.type.test = "logcounts", assay.type.ref = "logcounts")
saveRDS(bpe, paste0(outdir,prefix,"_bpe_main_annotated.RDS"))

bpe.fine <- SingleR(sce_norm, ref=bpe.ref, labels=bpe.ref$label.fine,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
saveRDS(bpe.fine, paste0(outdir,prefix,"_bpe_fine_main_annotated.RDS"))


# Inspect
length(unique(hpca$pruned.labels)) #use this for choosing optimal cluster number 
table(hpca$pruned.labels) 
length(unique(bpe$pruned.labels)) #use this for choosing optimal cluster number
table(bpe$pruned.labels) 
length(unique(bpe.fine$pruned.labels)) #use this for choosing optimal cluster number
table(bpe.fine$pruned.labels)

# Write table
write.csv(as.data.frame(table(hpca$pruned.labels)),paste0(outdir,prefix,"_cell_based_annotation_hpca_main.csv"))
write.csv(as.data.frame(table(bpe$pruned.labels)),paste0(outdir,prefix,"_cell_based_annotation_bpe_main.csv"))
write.csv(as.data.frame(table(bpe.fine$pruned.labels)),paste0(outdir,prefix,"_cell_based_annotation_bpe_fine.csv"))

# Update object
seurat_object@meta.data$hpca.main <- hpca$pruned.labels
seurat_object@meta.data$bpe.main <- bpe$pruned.labels
seurat_object@meta.data$bpe.fine <- bpe.fine$pruned.labels

# Transformed
if(processing.mode=="norm"){
  # Transform after annotation
  message("Normalizing and transforming object . . .")
  message("Processing normalized and transformed object . . .")
  seurat_object <- SCTransform(seurat_object, vst.flavor = "v2", verbose = FALSE)
}else{
  message("Processing normalized and transformed object . . .")
  seurat_object <- seurat_object
}

# Dimension Reduction
seurat_object <- RunPCA(seurat_object,features = rownames(seurat_object) , verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:50, verbose = FALSE)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:50, verbose = FALSE)
res.vec <- c(0.1,0.3,0.5,0.7,0.9)
for(res in res.vec){
  seurat_object <- FindClusters(seurat_object, resolution = res, verbose = FALSE)
}
seurat_object <- RunTSNE(object = seurat_object)

# Cell annotation plots
annos <- c("hpca.main","bpe.main","bpe.fine")
reductions <- c("umap","tsne")

#print
pdf(paste0(outdir,prefix,"_cell_based_annotation.pdf"))
for(anno in annos){
  for(red in reductions){
    if(length(unique(seurat_object[[anno]][[1]]))>20){
      plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = anno) + NoLegend()
    }else{
      plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = anno)
    }
    print(plot)
  }
}
dev.off()


# Cluster plots
clusters <- paste0(DefaultAssay(seurat_object),"_snn_res.",as.character(res.vec))
reductions <- c("umap","tsne")

# print
pdf(paste0(outdir,prefix,"_clusters.pdf"))
for(cluster in clusters){
  for(red in reductions){
    if(length(unique(seurat_object[[cluster]][[1]]))>20){
      plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = cluster) + NoLegend()
    }else{
      plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = cluster)
    }
    print(plot)
  }
}
dev.off()


# Factor plots
factor.groups <- c("donor","sample","treatment","timepoint")

# print
pdf(paste0(outdir,prefix,"_factor_groups.pdf"))
for(factor in factor.groups){
  for(red in reductions){
    if(length(unique(seurat_object[[factor]][[1]]))>20){
      plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = factor) + NoLegend()
    }else{
      plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = factor)
    }
    print(plot)
  }
}
dev.off()


# Cluster based annotaton

# hpca main
scores <- c()
for(cluster in clusters){
  scores <- c(scores, abs(length(unique(hpca$pruned.labels))-length(unique(seurat_object[[cluster]][[1]]))))
}
best_clust <- clusters[which.min(scores)]
anno.test <- SingleR(sce_norm, ref=hpca.ref, labels=hpca.ref$label.main, clusters=seurat_object[[best_clust]][[1]])
write.csv(as.data.frame(anno.test),paste0(outdir,prefix,"_",best_clust,"_based_annotation_hpca_main.csv"))

seurat_object[["hpca_clust_anno"]] <- "CELLS"
for(cl in unique(seurat_object[[best_clust]][[1]])){
  seurat_object$hpca_clust_anno[seurat_object[[best_clust]]==cl] <- anno.test[cl,"pruned.labels"]
}

pdf(paste0(outdir,prefix,"_hpca_clust_anno.pdf"))
for(red in reductions){
  if(length(unique(seurat_object[["hpca_clust_anno"]][[1]]))>20){
    plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = "hpca_clust_anno") + NoLegend()
  }else{
    plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = "hpca_clust_anno")
  }
  print(plot)
}  
dev.off()


# bpe main
scores <- c()
for(cluster in clusters){
  scores <- c(scores, abs(length(unique(bpe$pruned.labels))-length(unique(seurat_object[[cluster]][[1]]))))
}
best_clust <- clusters[which.min(scores)]
anno.test <- SingleR(sce_norm, ref=bpe.ref, labels=bpe.ref$label.main, clusters=seurat_object[[best_clust]][[1]])
write.csv(as.data.frame(anno.test),paste0(outdir,prefix,"_",best_clust,"_based_annotation_bpe_main.csv"))

seurat_object[["bpe_clust_anno"]] <- "CELLS"
for(cl in unique(seurat_object[[best_clust]][[1]])){
  seurat_object$bpe_clust_anno[seurat_object[[best_clust]]==cl] <- anno.test[cl,"pruned.labels"]
}

pdf(paste0(outdir,prefix,"_bpe_clust_anno.pdf"))
for(red in reductions){
  if(length(unique(seurat_object[["bpe_clust_anno"]][[1]]))>20){
    plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = "bpe_clust_anno") + NoLegend()
  }else{
    plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = "bpe_clust_anno")
  }
  print(plot)
}  
dev.off()

# bpe fine
scores <- c()
for(cluster in clusters){
  scores <- c(scores, abs(length(unique(bpe.fine$pruned.labels))-length(unique(seurat_object[[cluster]][[1]]))))
}
best_clust <- clusters[which.min(scores)]
anno.test <- SingleR(sce_norm, ref=bpe.ref, labels=bpe.ref$label.fine, clusters=seurat_object[[best_clust]][[1]])
write.csv(as.data.frame(anno.test),paste0(outdir,prefix,"_",best_clust,"_based_annotation_bpe_fine.csv"))

seurat_object[["bpe_fine_clust_anno"]] <- "CELLS"
for(cl in unique(seurat_object[[best_clust]][[1]])){
  seurat_object$bpe_fine_clust_anno[seurat_object[[best_clust]]==cl] <- anno.test[cl,"pruned.labels"]
}

pdf(paste0(outdir,prefix,"_bpe_fine_clust_anno.pdf"))
for(red in reductions){
  if(length(unique(seurat_object[["bpe_fine_clust_anno"]][[1]]))>20){
    plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = "bpe_fine_clust_anno") + NoLegend()
  }else{
    plot <- DimPlot(object = seurat_object, reduction = red, label = T, repel = T, group.by = "bpe_fine_clust_anno")
  }
  print(plot)
}  
dev.off()

# Save
saveRDS(seurat_object,paste0(outdir,prefix,"_processed_seurat.rds"))
#seurat_object <- readRDS(paste0(outdir,prefix,"_processed_seurat.rds"))





