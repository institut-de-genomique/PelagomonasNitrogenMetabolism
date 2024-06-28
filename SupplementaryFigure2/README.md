# Supplementary figure 2
This directory contains the files and script to generate Supplementary Figure 2 of the article.  
Input files are :  
Gene expression levels of *P. calceolata* RCC100 and 697 in different nitrate conditions : https://zenodo.org/uploads/12582059  

## Data loading
```{r}
TPM_table <- read.table("20230427_RCC100-Nitrate_transcriptomes_TPM.tsv")
```

## A) Cluster Dendrogram
```{r}
distance <- dist(t(TPM_table))
sup_fig2_A <- hclust(distance)
plot(sup_fig2_A)

```

## B) Correlation plot
```{r}
library(corrplot)
library(RColorBrewer)
corr_matrix <- cor(as.matrix(TPM_table))
cl <- c("#A52A2A","#FF0000","#E97451","#FAA0A0","#B9D9EB","#6699CC","#0066b2","#002D62")
col2 <- colorRampPalette(cl)
sup_fig2_B <- corrplot(corr_matrix, method="circle",is.corr = FALSE,col = col2(8) ,cl.lim = c(0,1)) 
sup_fig2_B
```
## Output
```{r}
#Data output
#write.table(as.data.frame(distance), file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/SupFigure2/output_supFig2_A",quote = F, sep="\t")
write.table(corr_matrix, file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/SupFigure2/output_supFig2_B",quote = F, sep="\t")

#Figure output
pdf("Sup_Figure_2_A.pdf", width = 9, height = 6)
print(plot(sup_fig2_A))
dev.off() 
pdf("Sup_Figure_2_B.pdf", width = 9, height = 6)
print(sup_fig2_B)
dev.off()
```
