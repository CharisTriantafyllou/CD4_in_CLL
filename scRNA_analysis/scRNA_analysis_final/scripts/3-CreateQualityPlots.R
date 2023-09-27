
################################# START ########################################

## Load libraries --------------------------------------------------------------

library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(dplyr)

## Load the Seurat Object ------------------------------------------------------

# Messages for OFF mode
if(filtering.mode=="OFF"){
  message("Filtering mode: OFF")
  prefix="preQC"
  message("Loading preQC object . . .")
  seurat_path <- paste0(obj.dir, prefix,"_seurat.rds")
  seurat_object <- readRDS(seurat_path)
  message("Creating PreQC plots . . .")
}

# Messages for ON mode
if(filtering.mode=="ON"){
  message("Filtering mode: ON")
  prefix="postQC"
  message("Loading postQC object . . .")
  seurat_path <- paste0(obj.dir, prefix,"_seurat.rds")
  seurat_object <- readRDS(seurat_path)
  message("Creating PostQC plots . . .")
}

## Make plots ------------------------------------------------------------------

# Define outdir
outdir=paste0(qual.plot.dir,prefix,"/")
dir.create(outdir,recursive = T)

# Mitochondrial genes depth percentage and gene number - create only for preQC object

if(filtering.mode=="OFF"){
  
  message("Creating mitochondrial percentage and feature number cut-off plots . . .")
  
  # Mitochondrial genes depth percentage 
  ggplot(seurat_object@meta.data, aes(y=percent.mt, x=nCount_RNA)) + geom_point() + geom_hline(yintercept=mito.cutoff.max)
  ggsave(paste0(outdir, prefix, "_MitoPlot.pdf"))
  
  # Number of genes across cells (many genes probably duplets)
  ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + geom_hline(yintercept=nGene.cutoff.max)
  ggsave(paste0(outdir, prefix, "_GenePlot.pdf"))
  
}else if(filtering.mode=="ON"){
  
  message("No need for mitochondrial percentage plots and feature number cut-off plots!")
  
}

# Violin Plots

v1 <- VlnPlot(seurat_object, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),
              ncol = 4,pt.size = 0) & theme(plot.title = element_text(size=10))
v2 <- VlnPlot(seurat_object, features = c("percent.IGH","percent.IGK","percent.IGL"), group.by = "timepoint",
              ncol = 3,pt.size = 0) & theme(plot.title = element_text(size=10))
v3 <- VlnPlot(seurat_object, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"), group.by = "timepoint",
        ncol = 4,pt.size = 0) & theme(plot.title = element_text(size=10))
v4 <- VlnPlot(seurat_object, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"), group.by = "treatment",
        ncol = 4,pt.size = 0) & theme(plot.title = element_text(size=10))

pdf(paste0(outdir, prefix,"_violin.pdf"))
print(v1)
print(v2)
print(v3)
print(v4)
dev.off()

# Scatterplots

fs1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "donor")
fs2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "donor")
fs3 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.rb", group.by = "donor")
fs4 <- FeatureScatter(seurat_object, feature1 = "percent.rb", feature2 = "percent.mt", group.by = "donor")
fs5 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample")
fs6 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
fs7 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.rb", group.by = "sample")
fs8 <- FeatureScatter(seurat_object, feature1 = "percent.rb", feature2 = "percent.mt", group.by = "sample")
fs9 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "timepoint")
fs10 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "timepoint")
fs11 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.rb", group.by = "timepoint")
fs12 <- FeatureScatter(seurat_object, feature1 = "percent.rb", feature2 = "percent.mt", group.by = "timepoint")
fs13 <- FeatureScatter(seurat_object, feature1 = "percent.rb", feature2 = "percent.mt", group.by = "Phase")

pdf(paste0(outdir, prefix,"_feature_correlation.pdf"))
print(fs1)
print(fs2)
print(fs3)
print(fs4)
print(fs5)
print(fs6)
print(fs7)
print(fs8)
print(fs9)
print(fs10)
print(fs11)
print(fs12)
print(fs13)
dev.off()

# Phase pie

ph1 <- as_tibble(seurat_object[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15)) + ggtitle(paste0(prefix," Cell Cycle scores"))

pdf(paste0(outdir, prefix,"_phase_pie.pdf"))
print(ph1)
dev.off()

# End

message("Done !!!")

################################## END #########################################



