
# Libraries
library(dplyr)
library(ggplot2)
library(enrichplot)
library(forcats)
library(ggstance)
library(EnhancedVolcano)

## ATAC-seq significant genes 

# Load 
TFs_up <- read.table("enrichR/Transcription_Factor_PPIs_ATAC_UP.txt",sep = "\t",header = T)
TFs_down <- read.table("enrichR/Transcription_Factor_PPIs_ATAC_DOWN.txt",sep = "\t",header = T)

# Extract gene number
TFs_up$Score <- sapply(TFs_up$Overlap,function(x)strsplit(x,split = "/")[[1]][1])
TFs_down$Score <- sapply(TFs_down$Overlap,function(x)strsplit(x,split = "/")[[1]][1])

# Edit
TFs_up$Score <- as.numeric(paste0("+",TFs_up$Score))
TFs_down$Score <- as.numeric(paste0("-",TFs_down$Score))

# Combine
tbl <- as_tibble(rbind(TFs_up,TFs_down))
write.csv(as.data.frame(tbl),"atac_enrichR.csv")

# Filter
tbl <- filter(tbl,Adjusted.P.value < 0.1)

# Re-arrange
y <- arrange(tbl, desc(Score)) %>% 
  group_by(sign(Score)) 

write.csv(as.data.frame(y),"atac_enrichR_sign.csv")

#plot
bar1 <- ggplot(y, aes(Score, fct_reorder(Term, Score), fill=Adjusted.P.value)) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL) + ggtitle("Changes in target accessibility per TF")

pdf("atac_enrichR.pdf")
bar1
dev.off()

write.csv(as.data.frame(y),"atac_enrichR.csv")

## ATAC-seq all results
library(EnhancedVolcano)
df <- read.csv("ATAC_CD4_NEW.csv",header = T,row.names = 1)

res <- as.data.frame(df[,c("SYMBOL","Fold","p.value","FDR")])

EnhancedVolcano(res,
                lab = res$SYMBOL,
                x = 'Fold',
                y = 'FDR',
                pCutoff = 0.05)

v1 <- EnhancedVolcano(res,
                lab = res$SYMBOL,
                x = 'Fold',
                y = 'FDR',
                title = 'ATAC-seq: CD4 with CLL vs CD4 isolated',
                subtitle = paste0("Motifs Up: ", sum(res$FDR < 0.05 & res$Fold > 0),
                                  ", Motifs Down: ", sum(res$FDR < 0.05 & res$Fold < 0)),
                pCutoff = 0.05,
                FCcutoff = 0.125,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                xlim = c(-2.5, 2.5),
                ylim = c(0,10))

pdf("atac_volcano.pdf")
v1
dev.off()
