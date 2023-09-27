
## Setup Directories -----------------------------------------------------------

# Working Directory
wd<-"/data/home/cntriantafyllou/CD4inCLL/Analysis/scRNA_analysis_final/"
dir.create(wd,recursive = T)
setwd(wd)

# Script Directory 
scr.dir <- paste0(wd,"scripts/")
dir.create(scr.dir,recursive = T)

## Provide universal parameters for filtering and QC plots ---------------------

# Parameters for filtering - Set these up even if filtering doesn't take place.
mito.cutoff.max <- 15 # keep cell with no more than 15% mitochondrial gene counts
nGene.cutoff.max <- 3000 # keep cell with no more than 3000 detected genes
nGene.cutoff.min <- 200 # keep cell with no less than 200 detected genes
min.cells.perc <- 0.05 # keep genes with non-zero count in at least 5% of cells
min.cells.n <- 3 # keep genes with non-zero count in at 3 of cells

## Choose modes and options ----------------------------------------------------

running.mode=7  # A number from 1 to 7
filtering.mode="ON" # ON or OFF, for running modes 1 and 3
cell_filt_type="number" # "number" or "percentage", for running mode 1
processing.mode="sct_integration" # For running mode 4, 5, 6, 7
use_integrated_assay=T # Only for 'sct_integration' mode

## Loading Data ----------------------------------------------------------------

if(running.mode==1){
  
  message("Mode 1: Creating Seurat Objects")
  
  # Filtering Mode
  #filtering.mode= "OFF"
  #filtering.mode= "ON"
  
  # Messages for OFF mode
  if(filtering.mode=="OFF"){
    message("Filtering mode: OFF")
    prefix="preQC"
  }
  
  # Messages for ON mode
  if(filtering.mode=="ON"){
    message("Filtering mode: ON")
    prefix="postQC"
  }
  
  # Object directory
  obj.dir <- paste0(wd,"seurat_objects/")
  dir.create(obj.dir,recursive = T)
  
  # Project Name
  project_name=paste0(prefix,"_Seurat")
  
  # Edit and/or Run Script #Deactivate after executing
  scr_name <- "1-CreateSeuratObjects.R"
  #file.create(paste0(scr.dir,scr_name))
  file.edit(paste0(scr.dir,scr_name))
  #source(paste0(scr.dir,scr_name))
  
}


## Batch effects ---------------------------------------------------------------

if(running.mode==2){
  
  message("Mode 2: Batch Effects Estimations")
  
  # Define object path
  obj.path <- paste0(obj.dir,"preQC","_seurat.rds")
  
  # Output directory
  batch.dir <- paste0(wd,"batch_effect_plots/")
  dir.create(batch.dir,recursive = T)
  
  # Edit and/or Run Script #Deactivate after executing
  scr_name <- "2-BatchEffects.R"
  #file.create(paste0(scr.dir,scr_name))
  file.edit(paste0(scr.dir,scr_name))
  #source(paste0(scr.dir,scr_name))
  
}
  

## Quality Plots -----------------------------------------------------------------

if(running.mode==3){
  
  message("Mode 3: Quality Plots")
  
  # Output directory
  qual.plot.dir=paste0(wd,"quality_plots/")
  dir.create(qual.plot.dir,recursive = T)
  
  # Edit and/or Run Script
  scr_name <- "3-CreateQualityPlots.R"
  #file.create(paste0(scr.dir,scr_name))
  file.edit(paste0(scr.dir,scr_name))
  #source(paste0(scr.dir,scr_name))
  
}


## PCA, Clustering, Annotation -------------------------------------------------

if(running.mode==4){
  
  message("Mode 4: PCA, Clustering and Annotation")
  
  # Processing mode
  #processing.mode="norm"
  #processing.mode="sct"
  #processing.mode="sct_sample_correction"
  #processing.mode="sct_all_corrections"
  #processing.mode="sct_integration"
  
  # Only for 'sct_integration' mode
  #use_integrated_assay=T
  #use_integrated_assay=F
  
  # Define object path
  obj.path <- paste0(obj.dir,"postQC","_seurat.rds")
  
  # Output directory
  tr_cl_ann.dir=paste0(wd,"trans_clust_anno/")
  dir.create(tr_cl_ann.dir,recursive = T)
  
  # Edit and/or Run Script
  scr_name <- "4-TransformationAndClustering.R"
  #file.create(paste0(scr.dir,scr_name))
  file.edit(paste0(scr.dir,scr_name))
  #source(paste0(scr.dir,scr_name))
  
}
  

## Annotation curation ---------------------------------------------------------

if(running.mode==5){
  
  message("Mode 5: Annotation curation")
  
  # Select the object that underwent this processing
  #processing.mode="sct_integration"
  #use_integrated_assay=T
  
  # Output/Working directory
  tr_cl_ann.dir=paste0(wd,"trans_clust_anno/")
  dir.create(tr_cl_ann.dir,recursive = T)
  
  # Object directory - where it is stored
  obj.dir=paste0(tr_cl_ann.dir,processing.mode,"/")
  
  # Edit and/or Run Script
  scr_name <- "5-AnnotationCuration.R"
  #file.create(paste0(scr.dir,scr_name))
  file.edit(paste0(scr.dir,scr_name))
  #source(paste0(scr.dir,scr_name))
  
  }


## TF activity -----------------------------------------------------------------

if(running.mode==6){
  
  message("Mode 6: TF Activity Analysis")
  
  # Select the object that underwent this processing
  #processing.mode="sct_integration"
  #use_integrated_assay=T
  
  # Object directory - where it is stored
  tr_cl_ann.dir=paste0(wd,"trans_clust_anno/")
  obj.dir=paste0(tr_cl_ann.dir,processing.mode,"/")
  
  # Output/Working directory
  outdir=paste0(wd,"decoupleR/")
  dir.create(outdir,recursive = T)
  
  # Edit and/or Run Script
  scr_name <- "6-decoupleR_TF.R"
  #file.create(paste0(scr.dir,scr_name))
  file.edit(paste0(scr.dir,scr_name))
  #source(paste0(scr.dir,scr_name))
  
  }


## Differential Expression -----------------------------------------------------

if(running.mode==7){
  
  message("Mode 7: Differential Expression")
  
  # Select the object that underwent this processing
  #processing.mode="sct_integration"
  #use_integrated_assay=T
  
  # Object directory - where it is stored
  tr_cl_ann.dir=paste0(wd,"trans_clust_anno/")
  obj.dir=paste0(tr_cl_ann.dir,processing.mode,"/")
  
  # Output/Working directory
  outdir=paste0(wd,"DEx/")
  dir.create(outdir,recursive = T)
  
  # Edit and/or Run Script
  scr_name <- "7-DEx.R"
  #file.create(paste0(scr.dir,scr_name))
  file.edit(paste0(scr.dir,scr_name))
  #source(paste0(scr.dir,scr_name))
  
  }
  



