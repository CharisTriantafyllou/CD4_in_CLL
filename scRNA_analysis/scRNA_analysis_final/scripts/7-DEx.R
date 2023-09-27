
## Load libraries --------------------------------------------------------------

# We load the required packages
library(Seurat)

# Library to create spredsheets
library(openxlsx)

# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
#library(pheatmap)
library(EnhancedVolcano)

# Messages for 'sct_integration' mode
if(processing.mode=="sct_integration"){
  
  message("Processing mode: integration after sct before annotation")
  prefix="sct_integration"
  
  if(use_integrated_assay==T){
    message("Using integrated assay . . .")
    prefix=paste0("integrated_",prefix)
  }else if(use_integrated_assay==F){
    message("Using SCT assay . . .")
    prefix=paste0("SCT_",prefix)
  }
  
  # Load Seurat object
  seurat_object <- readRDS(paste0(obj.dir,prefix,"_final_annotated_object.rds"))
  
}


#seurat_object <- subset(x = seurat_object, subset = final.anno %in% c("CD4+ T-cells"))
#seurat_object <- subset(x = seurat_object, subset = donor %in% c("CLL5"))


## Cells frequencies -----------------------------------------------------------

# Cells dataframes

data <- seurat_object

time.tab <- table(data$final.anno,data$timepoint)

donor.tab <- table(data$final.anno,data$donor)

treatment.tab <- table(data$final.anno,data$treatment)


# Save cell info 
setwd(outdir)

#timepoint
time.df <- as.data.frame(cbind(time.tab[,1],
                               time.tab[,2],
                               time.tab[,3],
                               time.tab[,4],
                               time.tab[,5]))

colnames(time.df) <- colnames(time.tab)
time.df <- cbind(rownames(time.df),time.df)
colnames(time.df)[1] <- "Cell Type"

time.df.perc <- as.data.frame(cbind(time.tab[,1]/sum(time.tab[,1])*100,
                                    time.tab[,2]/sum(time.tab[,2])*100,
                                    time.tab[,3]/sum(time.tab[,3])*100,
                                    time.tab[,4]/sum(time.tab[,4])*100,
                                    time.tab[,5]/sum(time.tab[,5])*100))

colnames(time.df.perc) <- colnames(time.tab)
time.df.perc <- cbind(rownames(time.df.perc),time.df.perc)
colnames(time.df.perc)[1] <- "Cell Type"

openxlsx::write.xlsx(time.df,"cells_and_time_df.xlsx")
openxlsx::write.xlsx(time.df.perc,"cells_and_time_df_perc.xlsx")


#donor
donor.df <- as.data.frame(cbind(donor.tab[,1],
                                donor.tab[,2],
                                donor.tab[,3],
                                donor.tab[,4]))

colnames(donor.df) <- colnames(donor.tab)
donor.df <- cbind(rownames(donor.df),donor.df)
colnames(donor.df)[1] <- "Cell Type"

donor.df.perc <- as.data.frame(cbind(donor.tab[,1]/sum(donor.tab[,1])*100,
                                     donor.tab[,2]/sum(donor.tab[,2])*100,
                                     donor.tab[,3]/sum(donor.tab[,3])*100,
                                     donor.tab[,4]/sum(donor.tab[,4])*100))

colnames(donor.df.perc) <- colnames(donor.tab)
donor.df.perc <- cbind(rownames(donor.df.perc),donor.df.perc)
colnames(donor.df.perc)[1] <- "Cell Type"

openxlsx::write.xlsx(donor.df,"cells_and_donor_df.xlsx")
openxlsx::write.xlsx(donor.df.perc,"cells_and_donor_df_perc.xlsx")

#treatment
treatment.df <- as.data.frame(cbind(treatment.tab[,1],
                                    treatment.tab[,2]))

colnames(treatment.df) <- colnames(treatment.tab)
treatment.df <- cbind(rownames(treatment.df),treatment.df)
colnames(treatment.df)[1] <- "Cell Type"

treatment.df.perc <- as.data.frame(cbind(treatment.tab[,1]/sum(treatment.tab[,1])*100,
                                         treatment.tab[,2]/sum(treatment.tab[,2])*100))

colnames(treatment.df.perc) <- colnames(treatment.tab)
treatment.df.perc <- cbind(rownames(treatment.df.perc),treatment.df.perc)
colnames(treatment.df.perc)[1] <- "Cell Type"

openxlsx::write.xlsx(treatment.df,"cells_and_treatment_df.xlsx")
openxlsx::write.xlsx(treatment.df.perc,"cells_and_treatment_df_perc.xlsx")


## Create labels and perform DEx -----------------------------------------------

data$final.anno.char <- gsub("[+]","",data$final.anno)
data$final.anno.char <- gsub("-","",data$final.anno.char)
data$final.anno.char <- gsub(" ","",data$final.anno.char)

data$Type_Treatment <- paste(data$final.anno.char,data$treatment,sep = "_")
Idents(data) <- "Type_Treatment"

# Prepare
data <- PrepSCTFindMarkers(data)

# DEx - Wilcoxon
no_treatment <- FindMarkers(data, assay = "SCT", ident.1 = "CD4Tcells_NO", ident.2 = "CD4Tcells_YES",logfc.threshold = 0)

# Re arrange
no_treatment$gene <- rownames(no_treatment)
no_treatment <- no_treatment[,c("gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj")]
no_treatment_0.05 <- no_treatment[no_treatment$p_val < 0.05,]

# View
head(no_treatment, n = 15)
head(no_treatment_0.05, n = 15)

# Save
write.xlsx(no_treatment,"DEx_CD4_treatment_NO_vs_YES.xlsx")
write.xlsx(no_treatment_0.05,"DEx_CD4_treatment_NO_vs_YES_filtered.xlsx")

Idents(data) <- "final.anno"
DefaultAssay(data) <- "SCT"

p1 <- FeaturePlot(data, features = c("HLA-DRA", "EEF1A1", "CD74"),
                  split.by = "treatment", max.cutoff = 3, cols = c("grey", "red"))

plots <- VlnPlot(data, features =  c("HLA-DRA", "EEF1A1", "CD74"), split.by = "treatment",
                 group.by = "final.anno", pt.size = 0, combine = FALSE)
p2 <- wrap_plots(plots = plots, ncol = 1)

pdf("top_hits.pdf")
print(p1)
print(p2)
dev.off()

#browseURL("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8516689/")
#browseURL("https://link.springer.com/article/10.1007/s11033-020-05736-5")


ench_volc <- EnhancedVolcano(no_treatment,
                             title = 'No treatment vs Ibrutinib Treatment in CD4+ Cells',
                             subtitle = paste0("L2FC: 0.01, p-val: 0.05"),
                             lab = rownames(no_treatment),
                             x = 'avg_log2FC',
                             y = 'p_val',
                             pCutoff = 0.05,
                             FCcutoff = 0.1,
                             pointSize = 3.0,
                             labSize = 6.0,
                             xlim = c(-1.2,1.2))

pdf("volcano.pdf", width = 7, height = 7)
print(ench_volc)
dev.off()

