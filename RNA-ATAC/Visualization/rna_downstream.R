
# Libraries
library(DESeq2)
library(openxlsx)
library(EnhancedVolcano)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(forcats)
library(ggstance)

## RNA-seq significant genes 

# Load genes 

dds <- readRDS("CD4_CLL_dds.rds")
res <- results(dds,alpha = 0.05)

# filter NAs
res <- res[rowSums(is.na(res))==0,]

## RNA-seq all results

res$gene <- rownames(res)
openxlsx::write.xlsx(res,"DEx_results.xlsx")

res_sign <- res[res$padj <= 0.05,]
openxlsx::write.xlsx(res_sign,"DEx_results_sign.xlsx")

res_sign_up <- res_sign[res_sign$log2FoldChange > 0,]
openxlsx::write.xlsx(res_sign_up,"DEx_results_sign_up.xlsx")

res_sign_down <- res_sign[res_sign$log2FoldChange < 0,]
openxlsx::write.xlsx(res_sign_down,"DEx_results_sign_down.xlsx")

## Plot
v1 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'RNA-seq: CD4 with CLL vs CD4 isolated',
                      subtitle = paste0("Genes Up: ", sum(res$padj < 0.05 & res$log2FoldChange > 0),
                                        " ( ",sum(res$padj < 0.05 & res$log2FoldChange > 0.5)," )",
                                        ", Genes Down: ", sum(res$padj < 0.05 & res$log2FoldChange < 0),
                                        " ( ", sum(res$padj < 0.05 & res$log2FoldChange < -0.5)," )"),
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 6.0,
                      col=c('black', 'black', 'black', 'red3'),
                      colAlpha = 1)

pdf("rna_volcano.pdf")
v1
dev.off()


## EnrichR


# Load 
TFs_up <- read.table("enrichR/Transcription_Factor_PPIs_RNA_UP.txt",sep = "\t",header = T)
TFs_down <- read.table("enrichR/Transcription_Factor_PPIs_RNA_DN.txt",sep = "\t",header = T)

# Extract gene number
TFs_up$Score <- sapply(TFs_up$Overlap,function(x)strsplit(x,split = "/")[[1]][1])
TFs_down$Score <- sapply(TFs_down$Overlap,function(x)strsplit(x,split = "/")[[1]][1])

# Edit
TFs_up$Score <- as.numeric(paste0("+",TFs_up$Score))
TFs_down$Score <- as.numeric(paste0("-",TFs_down$Score))

# Combine
tbl <- as_tibble(rbind(TFs_up,TFs_down))
write.csv(as.data.frame(tbl),"rna_enrichR.csv")

# Filter
tbl <- filter(tbl,Adjusted.P.value < 0.1)

# Re-arrange
y <- arrange(tbl, desc(Score)) %>% 
  group_by(sign(Score)) 

write.csv(as.data.frame(y),"rna_enrichR_sign.csv")

#plot
bar1 <- ggplot(y, aes(Score, fct_reorder(Term, Score), fill=Adjusted.P.value)) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL) + ggtitle("Changes in target expression per TF")

pdf("rna_enrichR.pdf")
bar1
dev.off()



