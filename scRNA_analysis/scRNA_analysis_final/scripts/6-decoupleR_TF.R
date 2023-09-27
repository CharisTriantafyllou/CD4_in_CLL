

## Load libraries --------------------------------------------------------------

## We load the required packages
library(Seurat)
library(decoupleR)

# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

# To save tables
library(openxlsx)


## Load Data -------------------------------------------------------------------

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

seurat_object <- subset(x = seurat_object, subset = final.anno %in% c("CD4+ T-cells"))
#seurat_object <- subset(x = seurat_object, subset = donor %in% c("CLL5"))

# Set wd
setwd(outdir)


## Normalization and Inspection ------------------------------------------------

# Edit
data <- seurat_object
Idents(data) <- "treatment"
DefaultAssay(data) <- "RNA"
data <- NormalizeData(data)

# Inspect
DimPlot(data, reduction = "umap", label = TRUE, group.by = "treatment", pt.size = 0.5) + NoLegend()
pdf("cd4_treatment_dimplot.pdf")
dp1 <- DimPlot(data, reduction = "umap", label = TRUE, group.by = "treatment", pt.size = 0.5) + NoLegend()
print(dp1)
dev.off()
#ggsave("cell_types_across_donors.pdf")

# Load net
net <- readRDS("net.RDS")
net

# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA@data)
#saveRDS(mat,"mat.RDS")
saveRDS(mat,"cd4_mat.RDS")
#saveRDS(mat,"cll5_cd4_mat.RDS")
#mat <- readRDS("mat.RDS")
mat <- readRDS("cd4_mat.RDS")
#mat <- readRDS("cll5_cd4_mat.RDS")


## Run UlM ---------------------------------------------------------------------

# Run ulm - TAKES A LOT
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)
acts

#saveRDS(acts,"acts.R")
saveRDS(acts,"cd4_acts.R")
#saveRDS(acts,"cll5_cd4_acts.R")

#acts <- readRDS("acts.R")
acts <- readRDS("cd4_acts.R")
#acts <- readRDS("cll5_cd4_acts.R")


# Extract ulm and store it in tfsulm in pbmc
data[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "tfsulm"

# Scale the data
data <- ScaleData(data)
data@assays$tfsulm@data <- data@assays$tfsulm@scale.data

#saveRDS(data,"data.tfsulm.RDS")
saveRDS(data,"cd4_data.tfsulm.RDS")
#saveRDS(data,"cll5_cd4_data.tfsulm.RDS")

#data <- readRDS("data.tfsulm.RDS")
data <- readRDS("cd4_data.tfsulm.RDS")
#data <- readRDS("cll5_cd4_data.tfsulm.RDS")


## Cluster Heatmap plots -------------------------------------------------------

## Do it manually to avoid surprises

#Idents(data) <- "donor"
#Idents(data) <- "sample"
#Idents(data) <- "timepoint"
Idents(data) <- "treatment"


{
  # Extract activities from object as a long dataframe
  df <- t(as.matrix(data@assays$tfsulm@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(data)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%
    summarise(mean = mean(score))
  
  # Get top tfs with more variable means across clusters
  n_tfs <- 50
  
  tfs <- df %>%
    group_by(source) %>%
    summarise(std = sd(mean)) %>%
    arrange(-abs(std)) %>%
    head(n_tfs) %>%
    pull(source)
  
  # Subset long data frame to top tfs and transform to wide matrix
  top_acts_mat <- df %>%
    filter(source %in% tfs) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('cluster') %>%
    as.matrix()
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))
  
  # Plot
  pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
}

  
#pdf("cd4_all_samples_donor_TF_clust.pdf")
#p <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
#print(p)
#dev.off()
  
#pdf("cd4_all_samples_samples_TF_clust.pdf")
#p <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
#print(p)
#dev.off()

#pdf("cd4_all_samples_timepoint_TF_clust.pdf")
#p <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
#print(p)
#dev.off()

pdf("cd4_all_samples_treatment_TF_clust.pdf")
p <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
print(p)
dev.off()
  
  
## Collect statistcs -----------------------------------------------------------
  
#find p-values
tfs_ts <- data@assays$tfsulm@data
idx.no <- data[["treatment"]]=="NO"
idx.yes <- data[["treatment"]]=="YES"
  
tfs <- rownames(tfs_ts)
  
ps <- rep(0,length(tfs))
for (i in 1:length(tfs)){
  w.t <- wilcox.test(tfs_ts[i,idx.no],tfs_ts[i,idx.yes])
  ps[i] <- w.t$p.value
}

adj.ps <- p.adjust(ps,"BH")
#tfs[adj.ps < 1e-100]

p_vals <- as.data.frame(cbind(ps,adj.ps))
diff_ps <- as.data.frame(cbind(tfs,p_vals))
rownames(diff_ps) <- tfs

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Make the difference df and merge with statistics
df_no <- as.data.frame(df[df$cluster=="NO",])
df_yes <- as.data.frame(df[df$cluster=="YES",])
rownames(df_no) <- df_no$source
rownames(df_yes) <- df_yes$source
df_yes <- df_yes[rownames(df_no),] #make sure the order is the same

df_combined <- as.data.frame(cbind(diff_ps[rownames(df_no),],df_yes[["mean"]],df_no[["mean"]]))
colnames(df_combined)[c(4,5)] <- c("mean.yes","mean.no")
df_combined$mean.diff <- df_combined$mean.no - df_combined$mean.yes
df_combined <- df_combined[,c(1,4,5,6,2,3)]
colnames(df_combined)[c(1,5,6)] <- c("TF","p.val","adj.p.val")

write.xlsx(df_combined,"cd4_all_samples_treatment_diff_TF_activities_res.xlsx")

# Choose top 100 or top 50
Idents(data) <- "treatment"

# Get top tfs with more variable means (greatest mean sd) across clusters
n_tfs <- 50 #or 100
  
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()
  
# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)

# Save
#df_combined_100 <- df_combined[tfs,]
df_combined_50 <- df_combined[tfs,]
#write.xlsx(df_combined_100,"cd4_treatment_top_100_sd_tf_differences.xlsx")
write.xlsx(df_combined_50,"cd4_treatment_top_50_sd_tf_differences.xlsx")
#plot_100_df <- as_tibble(df_combined_100)
plot_50_df <- as_tibble(df_combined_50)

# Keep for plots
plot_df <- df[df$source %in% tfs,]
colnames(plot_df)[3] <- "MeanScore"
#write.xlsx(df_combined_100,"cd4_treatment_top_100_sd_tf_activities.xlsx")
write.xlsx(df_combined_50,"cd4_treatment_top_50_sd_tf_activities.xlsx")
plot_df_no <- plot_df[plot_df$cluster=="NO",]
plot_df_yes <- plot_df[plot_df$cluster=="YES",]


## Plots -----------------------------------------------------------------------


# Plot
gp1 <- ggplot(plot_df_no, aes(x = reorder(source, MeanScore), y = MeanScore)) + 
  geom_bar(aes(fill = MeanScore), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        #axis.text.y = element_text(size =10, face= "bold"),
        axis.text.y = element_text(size =10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Activity score of TFs in CD4+ T cells before Treatment")


# Plot
gp2 <- ggplot(plot_df_yes, aes(x = reorder(source, MeanScore), y = MeanScore)) + 
  geom_bar(aes(fill = MeanScore), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Activity score of TFs in CD4+ T cells after Treatment")


gp3 <- ggplot(plot_df, aes(x = reorder(source, MeanScore), y = MeanScore)) + 
  geom_bar(aes(fill = MeanScore), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Activity score of TFs in CD4+ T cells")

# Plot TFs found significant in diffTF

diff_tf_res <- read.xlsx("CLLvsIsolated_Paired.results.summary.xlsx","Brief")
diff_tf_res$TF

plot_df_sign <- plot_df[plot_df$source %in% diff_tf_res$TF,]

gp4 <- ggplot(plot_df_sign, aes(x = reorder(source, MeanScore), y = MeanScore)) + 
  geom_bar(aes(fill = MeanScore), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Activity score of significant TFs in CD4+ T cells")


tfs.ord.lvls <- levels(reorder(plot_df$source, plot_df$MeanScore))
clrs <- rep("black",length(tfs.ord.lvls))
clrs[tfs.ord.lvls %in% diff_tf_res$TF] <- "red"

gp5 <- ggplot(plot_df, aes(x = reorder(source, MeanScore), y = MeanScore)) + 
  geom_bar(aes(fill = MeanScore), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10, colour = clrs),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold", colour = clrs),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("TF Pathways") + ggtitle("Activity score of significant TFs in CD4+ T cells")

#colnames(plot_100_df)[4] <- "MeanDiff"
colnames(plot_50_df)[4] <- "MeanDiff"

#gp6 <- ggplot(plot_100_df, aes(x = reorder(TF, MeanDiff), y = MeanDiff)) +  
gp6 <- ggplot(plot_50_df, aes(x = reorder(TF, MeanDiff), y = MeanDiff)) + 
  geom_bar(aes(fill = MeanDiff), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10, colour = clrs),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold", colour = clrs),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("TF Pathways") + ggtitle("Difference of TF Activity score in CD4+ T cells")


#pdf("cd4_100_max_sd_TFs_and_diffTF_TFs.pdf")
pdf("cd4_50_max_sd_TFs_and_diffTF_TFs.pdf", width = 12, height = 8)
print(gp1)
print(gp2)
print(gp3)
print(gp4)
print(gp5)
print(gp6)
dev.off()


## diff TF results

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% diff_tf_res$TF) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
ph1 <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)


# Save
plot_df <- df[df$source %in% diff_tf_res$TF,]
colnames(plot_df)[3] <- "score"
plot_df_no <- plot_df[plot_df$cluster=="NO",]
plot_df_yes <- plot_df[plot_df$cluster=="YES",]

#write.csv(plot_df,"cd4_all_samples_tf_activities.csv")
#write.csv(plot_df,"cll5_cd4_tf_activities.csv")
write.csv(plot_df,"diffTF_tfs_cd4_tf_activities.csv")

# Plot
gp1 <- ggplot(plot_df_no, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Activity score of diffTF TFs in CD4+ T cells before Treatment")


# Plot
gp2 <- ggplot(plot_df_yes, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Activity score of diffTF TFs in CD4+ T cells after Treatment")


gp3 <- ggplot(plot_df, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Activity score of diffTF TFs in CD4+ T cells")


df_combined_diffTF <- df_combined[df_combined$TF %in% diff_tf_res$TF,]
write.xlsx(df_combined_diffTF,"diffTF_tfs_cd4_tf_differences.csv")
plot_df_combined_diffTF <- as_tibble(df_combined_diffTF)
colnames(plot_df_combined_diffTF)[4] <- "MeanDiff"

gp4 <- ggplot(plot_df_combined_diffTF, aes(x = reorder(TF, MeanDiff), y = MeanDiff)) + 
  geom_bar(aes(fill = MeanDiff), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10),
        #axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Activity Difference of diffTF highlighted TFs in CD4+ T cells NO vs YES")


pdf("cd4_difftf_res_tf_activities.pdf", width = 9, height = 7)
print(ph1)
print(gp1)
print(gp2)
print(gp3)
print(gp4)
dev.off()
