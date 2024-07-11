# Figure 5
This directory contains the files and scripts to generate Figure 5 of the article : Heatmap with L2FC values on genes of interest.  
Input files are :  
- Table of gene names known to be involved in *Pelagomonas* nitrogen metabolism : gene_list.txt in this directory.
- Tables of DESeq2 results for RCC100 and RCC697 : [https://zenodo.org/records/12726053]

The following code is executed on R version 4.1.1  

## Required R libraries
```r
library(reshape2)
library(tidyverse)
```

## Load the data

##### Retrieve our genes of interest (our 53 genes involved in the nitrogen cycle)
```r
interesting_genes <- read.table("gene_list.txt", sep = "\t", header=T)
```

##### Retrieve Log2 Fold Changes (L2FC) and ajusted p-values (padj) from DESeq2 results
```r
##### RCC100
results_RCC100_400 <- read.table(file="DESeq-results_RCC100_Nitrate400", sep = "\t", header=T)[,c(2,5)]
results_RCC100_200 <- read.table(file="DESeq-results_RCC100_Nitrate200", sep = "\t", header=T)[,c(2,5)]
results_RCC100_Ammonium <- read.table(file="DESeq-results_RCC100_Ammonium", sep = "\t", header=T)[,c(2,5)]
results_RCC100_Cyanate <- read.table(file="DESeq-results_RCC100_Cyanate", sep = "\t", header=T)[,c(2,5)]
results_RCC100_Urea <- read.table(file="DESeq-results_RCC100_Urea", sep = "\t", header=T)[,c(2,5)]

##### RCC697
results_RCC697_200 <- read.table(file="DESeq-results_RCC697_Nitrate200", sep = "\t", header=T)[,c(2,5)]
results_RCC697_50 <- read.table(file="DESeq-results_RCC697_Nitrate50", sep = "\t", header=T)[,c(2,5)]

#### Environment
results_RCC100_tara <- read.table(file="DEG_nitrate_envi_complete.tab", sep = "\t", header=T)[,c(2,6)]
```

## Format the data
##### Add prefix with strain+condition to each column, in each dataframe
```r
colnames(results_RCC100_400) <- paste0("RCC100-400_", colnames(results_RCC100_400))
colnames(results_RCC100_200) <- paste0("RCC100-200_", colnames(results_RCC100_200))
colnames(results_RCC100_Ammonium) <- paste0("RCC100-Ammonium_", colnames(results_RCC100_Ammonium))
colnames(results_RCC100_Urea) <- paste0("RCC100-Urea_", colnames(results_RCC100_Urea))
colnames(results_RCC100_Cyanate) <- paste0("RCC100-Cyanate_", colnames(results_RCC100_Cyanate))
colnames(results_RCC697_200) <- paste0("RCC697-200_", colnames(results_RCC697_200))
colnames(results_RCC697_50) <- paste0("RCC697-50_", colnames(results_RCC697_50))
colnames(results_RCC100_tara) <- paste0("RCC100-env_", colnames(results_RCC100_tara))
```

##### Bind RCC100 results together, same for RCC697
```r
results_RCC100 <- cbind(results_RCC100_400,results_RCC100_200,results_RCC100_Ammonium,results_RCC100_Urea,results_RCC100_Cyanate)
results_RCC697 <- cbind(results_RCC697_200,results_RCC697_50)
#### RCC100 : 16518 rows, 10 columns  /  RCC697 : 14214 rows, 4 columns (3 genes are not in RCC100 results)
```

##### Merge RCC100 results with RCC697 and Tara results
```r
results_RCC100_RCC697 <- merge(results_RCC100, results_RCC697, by="row.names", all.x=T, all.y=T) # 16521 rows, 15 columns
results_RCC100_RCC697 <- merge(results_RCC100_RCC697, results_RCC100_tara, by.x="Row.names", by.y="row.names", all.x=T, all.y=T) # 16521 rows, 17 columns
```

##### Simplify colnames
```r
colnames(results_RCC100_RCC697) <- sub("log2FoldChange", "L2FC", colnames(results_RCC100_RCC697))
colnames(results_RCC100_RCC697)[colnames(results_RCC100_RCC697) == "Row.names"] <- "Gene"
```

## Merge with interesting genes table
```r
table0 <- merge(interesting_genes, results_RCC100_RCC697, by="Gene", all.x=F, all.y=F) # 52 rows, 17 columns
# Replace NA by 0
table0[is.na(table0)] <- 0
```

## Reshape heatmap data
```r
table0_melted <- table0 %>%
  pivot_longer(cols = -c("Gene", "Gene.name", "Function"),
               names_to = c("sample", ".value"), # only works with no majuscules
               names_sep = "_")
```
##### Add labels (to show asterisks)
```r
table0_melted <- table0_melted %>%
  mutate(label = case_when(
    (L2FC >= 2 | L2FC <= -2) & padj < 0.01 ~ "**",
    (L2FC >= 1 | L2FC <= -1) & padj < 0.01 ~ "*",
    padj >= 0.01 ~ "",
    TRUE ~ ""))
```
##### Cluster to order columns
```r
table0_melted[["sample"]] <- factor(table0_melted[["sample"]],levels=c("RCC100-env","RCC697-200","RCC697-50","RCC100-400","RCC100-200","RCC100-Ammonium","RCC100-Cyanate","RCC100-Urea"))
```

## Generate the heatmap
```r
# Heatmap with ggplot2
heatmap_tab1 <- ggplot(data = table0_melted, aes(x=sample, y=Gene, fill=L2FC, label=label)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid="white") +
  geom_text(nudge_y=-0.5, size=3) + # to show stars
  facet_grid(paste0(Function," (",Gene.name,")") ~., scales = "free", space = "free") +
  theme(strip.text.y  = element_text(angle=0, size=6), panel.background = element_blank(),
        panel.grid = element_blank(), axis.ticks = element_blank(), panel.spacing=unit(0.1, "lines"),
        axis.text.x = element_text(angle = 30, hjust=1, size=6),
        axis.text.y = element_text(size=5),
        legend.key.size = unit(0.2, 'cm'), 
        legend.key.height = unit(0.2, 'cm'), 
        legend.key.width = unit(0.2, 'cm'), 
        legend.title = element_text(size=7), 
        legend.text = element_text(size=6)) 
heatmap_tab1
```
