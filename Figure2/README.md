# Figure 2
This directory contains the files and scripts to generate Figure 2 of the article.  
Input files are :  
Gene expression levels of P. calceolata RCC100 and RCC697 in different nitrate conditions : https://zenodo.org/uploads/12582059  
Metadata associated to transcriptomes of RCC100 : metadata_transcriptomic_RCC100_nitrogen.tab in this directory  

The following code is executed on R version 4.1.1  

## Required R libraries
```{r}
library(DESeq2)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(tibble)
library(eulerr)
library(erer)
```

## Panels A to D
Differentially expressed genes of _P. calceolata_ RCC100 in 441 µM (A) and 220 µM (B) of nitrate and RCC697 in 220µM (C) and 50 µM (D) compared to 880 µM of nitrate. Genes with p-value < 0.01 and log2FC > 2 are coloured. 

### Input
```{r}
#for A and B
transcriptomic_count_RCC100_nitrogen <- read.table("20230427_RCC100-Nitrate_transcriptomes_rawcounts.tsv")
metadata <- read.table("metadata_transcriptomic_RCC100_nitrogen.tab")
#for C and D
results_RCC697_200 <- read.table("results_RCC697_200.tab")
results_RCC697_50 <- read.table("results_RCC697_50.tab")
```
### Figure 2A
```
#subset the count table
rawCountTable <- subset(transcriptomic_count_RCC100_nitrogen, select = c("NO3_800_R1", "NO3_800_R2", "NO3_800_R3"
                             , "NO3_400_R1", "NO3_400_R2", "NO3_400_R3"))
#subset the metadata
sampleInfo <- metadata %>% dplyr::filter(Condition %in% c("NO3_800", "NO3_400"))

#create DESeq2 object
se_star_matrix <- DESeqDataSetFromMatrix(countData = rawCountTable,
                                         colData = sampleInfo,
                                         design = ~ Condition)
#define control condition
se_star_matrix$Condition <- relevel(se_star_matrix$Condition, ref = "NO3_800")

#calculate DEG matrix
dds <- DESeq(se_star_matrix)

#shrink log2FC
results_RCC100_400 <- lfcShrink(dds, coef=resultsNames(dds)[2], type = 'normal')

#treat results to generate MA plot
figure_2A_output <- as.data.frame(results_RCC100_400)
figure_2A_output <- mutate(figure_2A_output,       
               regulation = case_when(
                 log2FoldChange > 2 & padj<0.01 ~ "downregulated"
                 ,         log2FoldChange < -2 & padj<0.01 ~ "upregulated"))

#extract top 10 most overexpressed and top 10 most underexpressed genes
top20_genes_2A <- figure_2A_output %>%
  rownames_to_column('Gene') %>%
  filter(regulation=="upregulated" | regulation=="downregulated") %>%
  group_by(regulation) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  dplyr::slice(1:10)
```

### Figure 2B
```{r}
#subset the count table
rawCountTable <- subset(transcriptomic_count_RCC100_nitrogen, select = c("NO3_800_R1", "NO3_800_R2", "NO3_800_R3"
                             , "NO3_200_R1", "NO3_200_R2", "NO3_200_R3"))

#subset the metadata
sampleInfo <- metadata %>% dplyr::filter(Condition %in% c("NO3_800", "NO3_200"))

#create DESeq2 object
se_star_matrix <- DESeqDataSetFromMatrix(countData = rawCountTable,
                                         colData = sampleInfo,
                                         design = ~ Condition)
#define control condition
se_star_matrix$Condition <- relevel(se_star_matrix$Condition, ref = "NO3_800")

#calculate DEG matrix
dds <- DESeq(se_star_matrix)

#shrink log2FC
results_RCC100_200 <- lfcShrink(dds, coef=resultsNames(dds)[2], type = 'normal')

#treat results to generate MA plot
figure_2B_output <- as.data.frame(results_RCC100_200)
figure_2B_output <- mutate(figure_2B_output,       
               regulation = case_when(
                 log2FoldChange > 2 & padj<0.01 ~ "downregulated"
                 ,         log2FoldChange < -2 & padj<0.01 ~ "upregulated"))

#extract top 10 most overexpressed and top 10 most underexpressed genes
top20_genes_2B <- figure_2B_output %>%
  rownames_to_column('Gene') %>%
  filter(regulation=="upregulated" | regulation=="downregulated") %>%
  group_by(regulation) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  dplyr::slice(1:10)
```
### Figure 2C
```{r}
figure_2C_output <- as.data.frame(results_RCC697_200)
figure_2C_output <- mutate(figure_2C_output, regulation = case_when(log2FoldChange > 2 & padj < 0.01 ~ "upregulated"
                                            ,log2FoldChange < -2 & padj < 0.01 ~ "downregulated"))
top20_genes_2C <- figure_2C_output %>%
  rownames_to_column('Gene') %>%
  filter(regulation=="upregulated" | regulation=="downregulated") %>%
  group_by(regulation) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  dplyr::slice(1:10)
top20_genes_2C

```

### Figure 2D  
```{r}
figure_2D_output <- as.data.frame(results_RCC697_50)
figure_2D_output <- mutate(figure_2D_output, regulation = case_when(log2FoldChange > 2 & padj < 0.01 ~ "upregulated"
                                            ,log2FoldChange < -2 & padj < 0.01 ~ "downregulated"))
top20_genes_2D <- figure_2D_output %>%
  rownames_to_column('Gene') %>%
  filter(regulation=="upregulated" | regulation=="downregulated") %>%
  group_by(regulation) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  dplyr::slice(1:10)
top20_genes_2D
```
## Panels E and F   
Euler diagram of DEGs in RCC100 and RCC697. The number of upregulated and downregulated genes is indicated.  

### Figure 2E  
```{r}
#subset results from DEG
up_DEG_400 <-subset(as.data.frame(results_RCC100_400),log2FoldChange > 2 & padj < 0.01)
down_DEG_400 <-subset(as.data.frame(results_RCC100_400),log2FoldChange < -2 & padj < 0.01)
up_DEG_200 <- subset(as.data.frame(results_RCC100_200), log2FoldChange > 2 & padj < 0.01)
down_DEG_200 <- subset(as.data.frame(results_RCC100_200), log2FoldChange < -2 & padj < 0.01)

#retrieve genes differentially expressed by conditions                  
figure_2E_output = 
  list(
    RCC100_800_vs_200_µM = c(rownames(up_DEG_200), rownames(down_DEG_200)),
    RCC100_800_vs_400_µM = c(rownames(up_DEG_400), rownames(down_DEG_400))
  )

####Number of gene for : 
##up-regulated in both condition
nrow(merge(up_DEG_400 , up_DEG_200, by="row.names", all=F))
##down-regulated in both condition
nrow(merge(down_DEG_400 , down_DEG_200, by="row.names", all=F))
#up-regulated in 200 µM only
length(setdiff(rownames(up_DEG_200), rownames(up_DEG_400)))
#down-regulated in 200 µM only
length(setdiff(rownames(down_DEG_200), rownames(down_DEG_400)))
#up-regulated in 400 µM only
length(setdiff(rownames(up_DEG_400), rownames(up_DEG_200)))
#down-regulated in 400 µM only
length(setdiff(rownames(down_DEG_400), rownames(down_DEG_200)))
```
### Figure 2F
```{r}
#subset results from DEG
up_DEG_50 <-subset(as.data.frame(results_RCC697_50),log2FoldChange > 2 & padj < 0.01)
down_DEG_50 <-subset(as.data.frame(results_RCC697_50),log2FoldChange < -2 & padj < 0.01)
up_DEG_200_697 <- subset(as.data.frame(results_RCC697_200), log2FoldChange > 2 & padj < 0.01)
down_DEG_200_697 <- subset(as.data.frame(results_RCC697_200), log2FoldChange < -2 & padj < 0.01)

#retrieve genes differentially expressed by conditions           
figure_2F_output = 
  list(
    RCC697_800_vs_200_µM = c(rownames(up_DEG_200_697), rownames(down_DEG_200_697)),
    RCC697_800_vs_50_µM = c(rownames(up_DEG_50), rownames(down_DEG_50))
  )

####Number of gene for : 
##up-regulated in both condition
nrow(merge(up_DEG_50 , up_DEG_200_697, by="row.names", all=F))
##down-regulated in both condition
nrow(merge(down_DEG_50 , down_DEG_200_697, by="row.names", all=F))
#up-regulated in 200 µM only
length(setdiff(rownames(up_DEG_200_697), rownames(up_DEG_50)))
#down-regulated in 200 µM only
length(setdiff(rownames(down_DEG_200_697), rownames(down_DEG_50)))
#up-regulated in 400 µM only
length(setdiff(rownames(up_DEG_50), rownames(up_DEG_200_697)))
#down-regulated in 400 µM only
length(setdiff(rownames(down_DEG_50), rownames(down_DEG_200_697)))
```

## Figure and Data output
```{r}
#Figure 2A
fig2A_MAplot <- ggplot()+  
  geom_point(data=figure_2A_output, aes(x=log2(baseMean), y=log2FoldChange, color=regulation), size=1)+  
  geom_text_repel(data=top20_genes_2A, aes(x=log2(baseMean), y=log2FoldChange, label=Gene), size=3,  point.size = NA) +
  ylim(-7.5,7.5) +  
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(-2, 2), color="red") +  
  labs(title = "DEG in changing nitrate concentration (882 µM vs 441 µM)", x="Mean of normalized count (TPM)", y="log2FC") +  
  theme_bw() +   
  scale_color_manual(values = c("#A3432A","#6B0F1F"), na.value = "grey80",name = "padj" , labels = c("overexpressed in 441 µM NO3","underexpressed in 441 µM NO3",  "non significative"))
fig2A_MAplot

#Figure_2B
fig2B_MAplot <- ggplot()+  
  geom_point(data=figure_2B_output, aes(x=log2(baseMean), y=log2FoldChange, color=regulation), size=1)+  
  geom_text_repel(data=top20_genes_2B, aes(x=log2(baseMean), y=log2FoldChange, label=Gene), size=3, point.size = NA) +
  ylim(-7.5,7.5) +  
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(-2, 2), color="red") +  
  labs(title = "DEG in changing nitrate concentration (882 µM vs 220 µM)", x="Mean of normalized count (TPM)", y="log2FC") +  
  theme_bw() +   
  scale_color_manual(values = c("#B86B37","#6B0F1F"), na.value = "grey80",name = "padj", labels = c("overexpressed in 220 µM NO3","underexpressed in 220 µM NO3",  "non significative"))
fig2B_MAplot

#Figure 2C
fig2C_MAplot <- ggplot()+  
  geom_point(data=figure_2C_output, aes(x=log2(baseMean), y=log2FoldChange, color=regulation), size=1)+  
  geom_text_repel(data=top20_genes_2C, aes(x=log2(baseMean), y=log2FoldChange, label=Gene), size=3, point.size = NA) + 
  ylim(-8.2,8.2) +  geom_hline(yintercept = 0) + 
  #geom_hline(yintercept = c(-2, 2), color="red") +  
  labs(title = "DEG in changing nitrate concentration (200 vs 800)", x="Normalized counts (median-of-ratios)", y="log2FC") +  
  theme_bw() +   
  scale_color_manual(values = c("#6B0F1F","#899D5F"), na.value = "grey80",name = "padj", labels = c("underexpressed in 200 µM","overexpressed in 200 µM",  "non significative"))
fig2C_MAplot

#Figure 2D
fig2D_MAplot <- ggplot()+  
  geom_point(data=figure_2D_output, aes(x=log2(baseMean), y=log2FoldChange, color=regulation), size=1, max.overlap=15)+  
  geom_text_repel(data=top20_genes_2D, aes(x=log2(baseMean), y=log2FoldChange, label=Gene), size=3, point.size = NA) + 
  ylim(-9.5,9.5) +  geom_hline(yintercept = 0) + 
  #geom_hline(yintercept = c(-2, 2), color="red") +  
  labs(title = "DEG in changing nitrate concentration (50 vs 800)", x="Normalized counts (median-of-ratios)", y="log2FC") +  
  theme_bw() +   
  scale_color_manual(values = c("#6B0F1F","#899D5F"), na.value = "grey80",name = "padj", labels = c("underexpressed in 50 µM","overexpressed in 50 µM",  "non significative"))
fig2D_MAplot

#Figure 2E
fig2E_euler <- plot(euler(figure_2E_output), quantities = TRUE, fills = c("#B86B37","#A3432A") #E5C038
     , main="RCC100 : DEG in low nitrate")
fig2E_euler

#Figure 2F
fig2F_euler <- plot(euler(figure_2F_output), quantities = TRUE, fills = c("#B86B37","#A3432A") #E5C038
     , main="RCC697 : DEG in low nitrate")
fig2F_euler

#Data output
write.table(figure_2A_output, file = "figure_2A_output.tab",quote = F, sep="\t")
write.table(figure_2B_output, file = "figure_2B_output.tab",quote = F, sep="\t")
write.table(figure_2C_output, file = "figure_2C_output.tab",quote = F, sep="\t")
write.table(figure_2D_output, file = "figure_2D_output.tab",quote = F, sep="\t")
write.list(figure_2E_output, file = "figure_2E_output.tab",quote = F)
write.list(figure_2F_output, file = "figure_2F_output.tab",quote = F)

#Figure output
#All figure are saved in one file
pdf("Figure_2.pdf", width = 9, height = 6)
print(fig2A_MAplot)
print(fig2B_MAplot)
print(fig2C_MAplot)
print(fig2D_MAplot)
print(fig2E_euler)
print(fig2F_euler)
dev.off() 
```

