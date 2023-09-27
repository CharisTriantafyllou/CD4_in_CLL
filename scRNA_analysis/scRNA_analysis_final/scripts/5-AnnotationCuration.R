
################################# START ########################################

## Load libraries --------------------------------------------------------------

library(Seurat) # Seurat related utilities
library(SeuratData) # Seurat related utilities
library(SeuratObject) # Seurat related utilities
library(sctransform) # SCTransform code and utilities
#library(glmGamPoi) # SCTransform model related code
#library(SingleCellExperiment) # SCE object related utilities
#library(SingleR) # Annotation package
#library(celldex) # Annotation related utilities
library(ggplot2) # Plot related utilities
library(dplyr) # Data wrangling
library(RColorBrewer) # For advanced color Palette options
library(patchwork) # For combining plots
#library(scater) # Optional single cell object utilities
#library(scuttle) # Optional single cell object utilities

# Messages for 'sct_integration' mode
if(processing.mode=="sct_integration"){
  
  message("Processing mode: integration after sct before annotation")
  prefix="sct_integration"
  
  # Define outdir
  outdir=paste0(tr_cl_ann.dir,prefix,"/")
  dir.create(outdir,recursive = T)
  
  if(use_integrated_assay==T){
    message("Using integrated assay . . .")
    prefix=paste0("integrated_",prefix)
  }else if(use_integrated_assay==F){
    message("Using SCT assay . . .")
    prefix=paste0("SCT_",prefix)
  }
  
  # Load Seurat object
  seurat_object <- readRDS(paste0(obj.dir,prefix,"_processed_seurat.rds"))
  
}

## Manual annotation ------------------------------------------------------------

# Some clusters like NK T cells, were pretty straighforward, but that was not always the case. 
# Especially differentiation between memory and naive T cells was impossible because of sparsity related issues.
# For this reason markers were estimated for each cluster and 

#marker_genes

for_plot <- seurat_object
DefaultAssay(for_plot) <- "SCT"

#for_markers <- seurat_object
#DefaultAssay(for_markers) <- "RNA"
#DE_markers <- FindAllMarkers(for_markers)
#saveRDS(DE_markers,paste0(obj.dir,prefix,"_0.9_cluster_markers.rds"))
DE_markers <- readRDS(paste0(obj.dir,prefix,"_0.9_cluster_markers.rds"))

DE_markers_6 <- DE_markers %>% 
  group_by(cluster) %>%
  slice_max(n = 6, order_by = avg_log2FC)

cl=2
FeaturePlot(for_plot, features = DE_markers_6$gene[DE_markers_6$cluster==cl])
VlnPlot(for_plot, features = DE_markers_6$gene[DE_markers_6$cluster==cl])

# Known markers from Hemberg's group scRNA-seq course
VlnPlot(for_plot,"IL7R") # CD4 T cells
VlnPlot(for_plot,"CD4") # CD4 T cells
VlnPlot(for_plot,"CD8A") # CD8 T cells
VlnPlot(for_plot,"CD8B") # CD8 T cells
VlnPlot(for_plot,"SELL") # Naive T cells
VlnPlot(for_plot,"CD44") # Memory T cells
VlnPlot(for_plot,"MS4A1") # B cells
VlnPlot(for_plot,"TCL1A") # CLL B cells - cancer related marker
#browseURL("https://www.proteinatlas.org/ENSG00000100721-TCL1A/single+cell+type") 

# Check this more inclusive list from the best practises book of Theis' lab
marker_genes = list(c("FCN1", "CD14"),
                    c("TCF7L2", "FCGR3A", "LYN"),
                    c("CD14","ID2","VCAN", "S100A9","CLEC12A","KLF4", "PLAUR"),
                    c("CLEC9A", "CADM1"),
                    c("CST3","COTL1","LYZ","DMXL2","CLEC10A","FCER1A"), # Note: DMXL2 should be negative
                    c("SLC4A1", "SLC25A37", "HBB", "HBA2", "HBA1", "TFRC"),
                    c("MKI67", "HBA1", "HBB"),
                    c("CDK6","SYNGR1","HBM","GYPA"), # Note HBM and GYPA are negative markers
                    c("GNLY", "NKG7", "CD247", "GRIK4", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"),
                    c("ID2", "PLCG2", "GNLY", "SYNE1"),
                    c("VPREB1","MME","EBF1","SSBP2","BACH2","CD79B","IGHM","PAX5","PRKCE","DNTT","IGLL1"),
                    c("MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"),
                    c("MS4A1","SSPN","ITGB1","EPHA4","COL4A4","PRDM1","IRF4","CD38",
                      "XBP1","PAX5","BCL11A","BLK","IGHD","IGHM","ZNF215"), # Note IGHD and IGHM are negative markers
                    c("MME", "CD38", "CD24", "ACSM3", "MSI2"),
                    c("MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"),
                    c("XBP1", "RF4", "PRDM1", "PAX5"), # Note PAX5 is a negative marker
                    c("CD4", "IL7R", "TRBC2", "ITGB1"),
                    c("CD4", "IL7R", "TRBC2", "CCR7"),
                    c("CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"),
                    c("CD69", "CD38"),  # CD69 much better marker!
                    c("LEF1", "CCR7", "TCF7"),
                    c("GZMB", "IL3RA", "COBLL1", "TCF4"),
                    c("MPO", "BCL2", "KCNQ5", "CSF3R"),
                    c("NRIP1", "MECOM", "PROM1", "NKAIN2", "CD34"),
                    c("ZNF385D","ITGA2B","RYR3","PLCB1")  # Note PLCB1 is a negative marker
) 

names(marker_genes) <- c("CD14+ Mono","CD16+ Mono","ID2-hi myeloid prog","cDC1","cDC2",
                         "Normoblast","Erythroblast", "Proerythroblast","NK","ILC","Lymph prog",
                         "Naive CD20+ B","B1 B","Transitional B","Plasma cells","Plasmablast",
                         "CD4+ T activated","CD4+ T naive","CD8+ T","T activation","T naive",
                         "pDC","G/M prog","HSC","MK/E prog") 


# Check each category separately because code didn't work as expected

#pdf(paste0(outdir,prefix,"_celltype_markers.pdf"))
#for(celltype in names(marker_genes)){
  celltype=names(marker_genes)[14]
  print(celltype)
  VlnPlot(for_plot, features = marker_genes[[celltype]])
  #FeaturePlot(for_plot, features = marker_genes[[celltype]])
  #FeaturePlot(seurat_object, features = marker_genes[[celltype]], assay)
#}
#dev.off()

# Visualization of clusters and previous annotation
DimPlot(for_plot,label=T ) #Identities are integrated_snn_res.0.9
#DimPlot(for_plot,label=T, group.by = "integrated_snn_res.0.9")
#DimPlot(for_plot,label=T, reduction ="tsne", group.by = "integrated_snn_res.0.9")

DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "bpe_fine_clust_anno") #Best annotation
DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "bpe_clust_anno") #No CD4 Tcells
DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "hpca_clust_anno") #No discrimination of Tcells categories
  
seurat_object$final.anno <- "cell"
seurat_object$final.anno[Idents(seurat_object) %in% c(5,13,20)] <- "Monocytes" 
seurat_object$final.anno[Idents(seurat_object) %in% c(12)] <- "NK T-cells" 
seurat_object$final.anno[Idents(seurat_object) %in% c(3,8,10,14)] <- "CD8+ T-cells"
seurat_object$final.anno[Idents(seurat_object) %in% c(11)] <- "CD4+ T-cells"
seurat_object$final.anno[Idents(seurat_object) %in% c(0,1,2,4,5,6,7,9,15,16,17,18,19)] <- "CLL B-cells"

# Save now
saveRDS(seurat_object,paste0(outdir,prefix,"_final_annotated_object.rds"))
write.csv(as.data.frame(table(seurat_object$final.anno)),paste0(outdir,prefix,"_final_cell_anno.csv"))

# Check also these based on the Rendeiro et al. paper:
#browseURL("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-14081-6/MediaObjects/41467_2019_14081_MOESM1_ESM.pdf")

#B cells (CD19+CD5-CD38-)
#CD4+ T cells (CD3+CD4+)
#CD8+ T cells (CD3+CD8+)
#CLL cells (CD19+CD5+),
#gamma-delta T-cells (CD3+CD4-CD8-),
#myeloid cells (CD3-CD19-CD14+),
#NK cells (CD3-CD19-CD56+),
#NK-T cells (CD3+CD8+CD56+),
#and T cells (CD3+) from PBMCs of patients 

FeaturePlot(for_plot, features = c("CD3D","CD8A","CD8B","NCAM1"))
VlnPlot(for_plot, features = c("CD3D","CD8A","CD8B","NCAM1"))

pdf(paste0(outdir,prefix,"_final_annotations.pdf"))

DimPlot(seurat_object, label=T, repel=T, reduction="umap" )
DimPlot(seurat_object, label=T, repel=T, reduction="tsne" )

DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "final.anno") + ggtitle("Final annotation UMAP")
DimPlot(seurat_object, label=T, repel=T, reduction="tsne", group.by = "final.anno") + ggtitle("Final annotation t-SNE")

DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "bpe.fine") + NoLegend() #mixed signals - too many labels
DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "bpe.main") + NoLegend()#mixed signals - too many labels
DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "hpca.main") + NoLegend()  #mixed signals - too many labels

DimPlot(seurat_object, label=T, repel=T, reduction="tsne", group.by = "bpe.fine") + NoLegend()#mixed signals - too many labels
DimPlot(seurat_object, label=T, repel=T, reduction="tsne", group.by = "bpe.main") + NoLegend()#mixed signals - too many labels
DimPlot(seurat_object, label=T, repel=T, reduction="tsne", group.by = "hpca.main") + NoLegend() #mixed signals - too many labels

DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "bpe_fine_clust_anno") #Best annotation
DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "bpe_clust_anno") #No CD4 Tcells
DimPlot(seurat_object, label=T, repel=T, reduction="umap", group.by = "hpca_clust_anno") #No discrimination of Tcells categories

DimPlot(seurat_object, label=T, repel =T, reduction="tsne", group.by = "bpe_fine_clust_anno") #Best annotation
DimPlot(seurat_object, label=T, repel =T, reduction="tsne", group.by = "bpe_clust_anno") #No CD4 Tcells
DimPlot(seurat_object, label=T, repel =T, reduction="tsne", group.by = "hpca_clust_anno") #No discrimination of Tcells categories

dev.off()

table(seurat_object$final.anno)


################################# START ########################################
