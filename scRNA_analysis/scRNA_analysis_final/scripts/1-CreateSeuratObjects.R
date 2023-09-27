
################################# START ########################################

## Load libraries --------------------------------------------------------------

library(Seurat)
library(SeuratData)
library(SeuratObject)

## Load the Expression Matrix --------------------------------------------------

# Load data and create the matrix
message("Creating the expression matrix . . . ")
data_dir=paste0(wd,"data/GSE111014") # Define data path 
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir) # Load data retrieved from GEO

## Create Seurat Object and Inspect --------------------------------------------

if(filtering.mode=="OFF"){
  
  message("Creating unfiltered seurat object . . .")
  seurat_object = CreateSeuratObject(counts = expression_matrix,
                                     min.cells = 0,
                                     min.features = 0,
                                     project = project_name)
  
  }else if(filtering.mode=="ON"){
    
    message("Creating filtered seurat object . . .")
     
    if(cell_filt_type == "number"){
      
      seurat_object = CreateSeuratObject(counts = expression_matrix,
                                         min.cells = min.cells.n,
                                         min.features = nGene.cutoff.min,
                                         project = project_name)
      
    }else if(cell_filt_type == "percentage"){
      
      seurat_object = CreateSeuratObject(counts = expression_matrix,
                                         min.cells = round(min.cells.perc*ncol(expression_matrix)),
                                         min.features = 0,
                                         project = project_name)

    }
    
    }

# Inspect object 
message(paste0("Dataset consists of ", nrow(seurat_object)," genes and ", ncol(seurat_object), " cells"))
sce <- as.SingleCellExperiment(seurat_object)
tmpcount <- counts(sce)
avg.depth=mean(tmpcount[tmpcount>0])
colsum = colSums(tmpcount == 0)
sparsity = sum(colsum)/(nrow(sce)*ncol(sce))
metrics = c(max(tmpcount), round(avg.depth, 2), round(sparsity, 2))
metric.names = c("Max count: ", "Average Depth (count): ","Sparsity: ")
message(paste0("Max count: ", max(tmpcount)))
message(paste0("Average Depth (count): ", round(avg.depth, 2)))
message(paste0("Sparsity: ", round(sparsity, 2)))
csv <- as.data.frame(cbind(metric.names, metrics))
write.csv(csv, paste0(obj.dir,prefix,"_depth_metrics.csv"))
rm(sce) # remove sce object

## Extract information data from cell names ------------------------------------

meta <- seurat_object@meta.data

# Define functions
barcodes <- rownames(meta)
find_donor <- function(x){
  donor <- strsplit(strsplit(x,"-")[[1]][3],"_")[[1]][2]
  return(donor)
}
find_sample <- function(x){
  sample <- paste(strsplit(x,"_")[[1]][2:3],collapse="_")
  return(sample)
}
find_condition <- function(x){
  condition <- strsplit(strsplit(x,"-")[[1]][3],"_")[[1]][3]
  return(condition)
}

# Apply functions
donors <- sapply(barcodes, find_donor)
samples <- sapply(barcodes, find_sample)
conditions <- sapply(barcodes, find_condition)
days <- gsub("^d","",conditions)
barcode_df <- as.data.frame(cbind(barcodes,donors,samples,days))
rownames(barcode_df) <- NULL
barcode_df$days <- as.numeric(barcode_df$days)
unique(barcode_df$donors)

## Add experiment information to metadata --------------------------------------

seurat_object[["orig.ident"]] <- "ibrutinib_scRNAseq"
seurat_object[["sample"]] <- barcode_df$samples
seurat_object[["donor"]] <- barcode_df$donors
seurat_object[["timepoint"]] <- barcode_df$days
seurat_object[["treatment"]] <- "YES"
seurat_object[["treatment"]][seurat_object[["timepoint"]]==0] <- "NO"
sce <- as.SingleCellExperiment(seurat_object)

## Add gene percentages information --------------------------------------------

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-") #percentage of counts corresponding to mitochondrial genes
seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")

## Second round of Filtering  --------------------------------------------------

if(filtering.mode=="ON"){
  
  message("Second round of filtering . . .")
  
  # IG genes
  allGenes <- rownames(seurat_object)
  ighGenes <- allGenes[(grepl("^IGH", allGenes) & !grepl("^IGHMBP", allGenes)) | grepl("^IGL", allGenes) | grepl("^IGK", allGenes)]
  seurat_object <- seurat_object[allGenes[!allGenes %in% ighGenes],]
  
  message("IGH gemes removed!")
  
  # Gene number and Mitochondrial number 
  seurat_object <- subset(seurat_object, subset = nFeature_RNA <= nGene.cutoff.max)
  seurat_object <- subset(seurat_object, subset = percent.mt <= mito.cutoff.max)
  
  message("Cells with high gene content and mitochondrial gene content removed!")
  
}else if(filtering.mode=="OFF"){
  message("No second round of filtering performed")
}


## Add gene percentages information --------------------------------------------

seurat_object[["percent.IGH"]] <- PercentageFeatureSet(seurat_object, pattern = "^IGH")
seurat_object[["percent.IGK"]] <- PercentageFeatureSet(seurat_object, pattern = "^IGK")
seurat_object[["percent.IGL"]] <- PercentageFeatureSet(seurat_object, pattern = "^IGL")

## Perform Cell cycle scoring --------------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
for_cc_scoring <- NormalizeData(seurat_object) # only count normalization is needed here 
cc_scored <- CellCycleScoring(for_cc_scoring,
                              s.features = s.genes,
                              g2m.features = g2m.genes)
 
## Add cell cycle information by updating metadata -----------------------------

seurat_object@meta.data <- cc_scored@meta.data # update metadata
rm(cc_scored) # remove object

## Save object -----------------------------------------------------------------

message("Saving seurat object to: ", paste0(obj.dir, prefix,"_seurat.rds"))
saveRDS(seurat_object, paste0(obj.dir,prefix,"_seurat.rds"))

################################## END #########################################
