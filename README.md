---
title: "SupFigure2"
output: html_document
date: '2024-06-13'
---

Supplementary figure 1 

Data preparation
```{r}
TPM_table <- read.table("/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/SupFigure2/transcriptomic_TPM_RCC100_nitrogen.tab")
```

A) Cluster Dendrogram
```{r}
distance <- dist(t(TPM_table))
sup_fig2_A <- hclust(distance)
plot(sup_fig2_A)

```

B) Correlation plot
```{r}
library(corrplot)
library(RColorBrewer)
corr_matrix <- cor(as.matrix(TPM_table))
cl <- c("#A52A2A","#FF0000","#E97451","#FAA0A0","#B9D9EB","#6699CC","#0066b2","#002D62")
col2 <- colorRampPalette(cl)
sup_fig2_B <- corrplot(corr_matrix, method="circle",is.corr = FALSE,col = col2(8) ,cl.lim = c(0,1)) 
sup_fig2_B
```
Output
```{r}
#Data output
#write.table(as.data.frame(distance), file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/SupFigure2/output_supFig2_A",quote = F, sep="\t")
write.table(corr_matrix, file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/SupFigure2/output_supFig2_B",quote = F, sep="\t")

#Figure output
#All figure are saved in one file
pdf("/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/SupFigure2/Sup_Figure_2_A.pdf", width = 9, height = 6)
print(plot(sup_fig2_A))
dev.off() 

pdf("/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/SupFigure2/Sup_Figure_2_B.pdf", width = 9, height = 6)
print(sup_fig2_B)
dev.off()
```
