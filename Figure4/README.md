# Figure 4
Figure 4 shows the transcriptomic response of P. calceolata RCC100 cultivated with different nitrogen compounds. Following codes generates MA-plot and Venn diagram from transciptomic table for RCC100.

20230427_RCC100-Nitrate_transcriptomes_rawcounts.tsv is available here : https://zenodo.org/records/12582059  
metadata_transcriptomic_RCC100-RCC697_nitrogen.tab is [here](https://github.com/institut-de-genomique/PelagomonasNitrogenMetabolism/blob/main/Figure2/metadata_transcriptomic_RCC100-RCC697_nitrogen.tsv).  

A, B, C) Differentially expressed genes in 882 µM of ammonium (B), 441 µM of urea (C) and 882 µM of cyanate (D) compared to 882 µM of nitrate. Genes with p-value < 0.01 and log2FC > 2 are coloured. 
D left, Dright) Euler diagram of genes overexpressed (E) or under expressed (F) in at least one of the alternative nitrogen sources.  

## Library loading
```r
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tibble)
library(eulerr)
library(erer)
```

## Input
```r
transcriptomic_count_RCC100_nitrogen <- read.table("20230427_RCC100-Nitrate_transcriptomes_rawcounts.tsv")

metadata <- read.table("metadata_transcriptomic_RCC100-RCC697_nitrogen.tab")
```

## Treatment
```r
#subset the count table
rawCountTable <- subset(transcriptomic_count_RCC100_nitrogen, select = c(
                              "NO3_800_R1", "NO3_800_R2", "NO3_800_R3"
                             , "NH4_R1", "NH4_R2", "NH4_R3"
                             , "Urea_R1", "Urea_R2", "Urea_R3"
                             , "Cyanate_R1", "Cyanate_R2", "Cyanate_R3"))
#subset the metadata
sampleInfo <- metadata %>% dplyr::filter(Condition %in% c("NO3_800", "NH4", "Urea", "Cyanate"))


#create DESeq2 object
se_star_matrix <- DESeqDataSetFromMatrix(countData = rawCountTable,
                                         colData = sampleInfo,
                                         design = ~ Source)

#define control condition
se_star_matrix$Source <- relevel(se_star_matrix$Source, ref = "NO3")

#calculate DEG matrix
dds <- DESeq(se_star_matrix)

#shrink log2FC
res2_NO3_Cyanate <- lfcShrink(dds, coef="Source_Cyanate_vs_NO3", type="apeglm")
res2_NO3_Urea <- lfcShrink(dds, coef="Source_Urea_vs_NO3", type="apeglm")
res2_NO3_NH4 <- lfcShrink(dds, coef="Source_NH4_vs_NO3", type="apeglm")

```

Figure 4A
```r
#treat result to generate MA plot
figure_4A_output <- as.data.frame(res2_NO3_NH4)
figure_4A_output <- mutate(figure_3B_output, regulation = case_when(log2FoldChange > 2 & padj<0.01 ~ "upregulated"
                                            ,log2FoldChange < -2 & padj<0.01 ~ "downregulated"))

#extract top 10 most overexpressed and top 10 most underexpressed genes
top20_genes_3B <- figure_3B_output %>%
  rownames_to_column('Gene') %>%
  filter(regulation=="upregulated" | regulation=="downregulated") %>%
  group_by(regulation) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  dplyr::slice(1:10)
```


Figure 4B
```r
#treat result to generate MA plot
figure_4B_output <- as.data.frame(res2_NO3_Urea)
figure_4B_output <- mutate(figure_3C_output, regulation = case_when(log2FoldChange > 2 & padj<0.01 ~ "upregulated"
                                            ,log2FoldChange < -2 & padj<0.01 ~ "downregulated"))

#extract top 10 most overexpressed and top 10 most underexpressed genes
top20_genes_3C <- figure_3C_output %>%
  rownames_to_column('Gene') %>%
  filter(regulation=="upregulated" | regulation=="downregulated") %>%
  group_by(regulation) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  dplyr::slice(1:10)

```

Figure 4C
```r
#treat result to generate MA plot
figure_4C_output <- as.data.frame(res2_NO3_Cyanate)
figure_4C_output <- mutate(figure_3D_output, regulation = case_when(log2FoldChange > 2 & padj<0.01 ~ "upregulated"
                                            ,log2FoldChange < -2 & padj<0.01 ~ "downregulated"))

#extract top 10 most overexpressed and top 10 most underexpressed genes
top20_genes_3D <- figure_3D_output %>%
  rownames_to_column('Gene') %>%
  filter(regulation=="upregulated" | regulation=="downregulated") %>%
  group_by(regulation) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  dplyr::slice(1:10)

```


Figure 4D left
```r
up_DEG_Cyanate <- subset(as.data.frame(res2_NO3_Cyanate),log2FoldChange > 2 & padj < 0.01)
up_DEG_NH4 <- subset(as.data.frame(res2_NO3_NH4),log2FoldChange > 2 & padj < 0.01)
up_DEG_Urea <- subset(as.data.frame(res2_NO3_Urea),log2FoldChange > 2 & padj < 0.01)

figure_4Dl_output = 
  list(
    Cyanate=rownames(up_DEG_Cyanate),
    NH4=rownames(up_DEG_NH4),
    Urea=rownames(up_DEG_Urea)
  )
```

Figure 4D right
```r
down_DEG_Cyanate <- subset(as.data.frame(res2_NO3_Cyanate),log2FoldChange < -2 & padj < 0.01)
down_DEG_NH4 <- subset(as.data.frame(res2_NO3_NH4),log2FoldChange < -2 & padj < 0.01)
down_DEG_Urea <- subset(as.data.frame(res2_NO3_Urea),log2FoldChange < -2 & padj < 0.01)

figure_4Dr_output = 
  list(
    Cyanate=rownames(down_DEG_Cyanate),
    NH4=rownames(down_DEG_NH4),
    Urea=rownames(down_DEG_Urea)
  )
```

## Output
```r
#figure 4A
fig4A_MAplot <- ggplot()+  
  geom_point(data=figure_4A_output, aes(x=log2(baseMean), y=log2FoldChange, color=regulation), size=1)+  
  geom_text_repel(data=top20_genes_4A, aes(x=log2(baseMean), y=log2FoldChange, label=Gene), size=3, point.size = NA) + 
  ylim(-7.5,7.5) +  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(-2, 2), color="red") +  
  labs(title = "DEG in changing nitrogen source (NO3 vs ammonium)", x="Mean of normalized count (TPM)", y="log2FC") +  
  theme_bw() +   
  scale_color_manual(values = c("#6B0F1F", "#E5C038"), na.value = "grey80",name = "padj", labels = c("underexpressed in ammonium","overexpressed in ammonium",  "non significative"))
fig4A_MAplot

#figure 4B
fig4B_MAplot <- ggplot()+  
  geom_point(data=figure_4B_output, aes(x=log2(baseMean), y=log2FoldChange, color=regulation), size=1)+  
  geom_text_repel(data=top20_genes_4B, aes(x=log2(baseMean), y=log2FoldChange, label=Gene), size=3, point.size = NA) + 
  ylim(-8.2,8.2) +  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(-2, 2), color="red") +  
  labs(title = "DEG in changing nitrogen source (NO3 vs urea)", x="Mean of normalized count (TPM)", y="log2FC") +  
  theme_bw() +   
  scale_color_manual(values = c("#6B0F1F","#899D5F"), na.value = "grey80",name = "padj", labels = c("underexpressed in urea","overexpressed in urea",  "non significative"))
fig4B_MAplot

#figure 4C
fig3D_MAplot <- ggplot()+  
  geom_point(data=figure_4C_output, aes(x=log2(baseMean), y=log2FoldChange, color=regulation), size=1)+  
  geom_text_repel(data=top20_genes_4C, aes(x=log2(baseMean), y=log2FoldChange, label=Gene), size=3, point.size = NA) + 
  ylim(-10,6) +  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(-2, 2), color="red") +  
  labs(title = "DEG in changing nitrogen source (NO3 vs cyanate)", x="Mean of normalized count (TPM)", y="log2FC") +  
  theme_bw() +   
  scale_color_manual(values = c("#6B0F1F","#297373"), na.value = "grey80",name = "padj", labels = c("underexpressed in cyanate","overexpressed in cyanate",  "non significative"))
fig4C_MAplot

#figure 4D left
fig4Dl_euler <- plot(euler(figure_4Dl_output), quantities = TRUE, fills = c("#297373","#E5C038", "#899D5F")
     , main="Underexpressed in alternative nitrogen sources")
fig4Dl_euler

#figure 4D right
fig4Dr_euler <- plot(euler(figure_4Dr_output), quantities = TRUE, fills = c("#297373","#E5C038", "#899D5F")
     , main="Underexpressed in alternative nitrogen sources")
fig4Dr_euler


#Data output
write.table(figure_4A_output, file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/Figure3/figure_3B_output.tab",quote = F, sep="\t")
write.table(figure_4B_output, file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/Figure3/figure_3C_output.tab",quote = F, sep="\t")
write.table(figure_4C_output, file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/Figure3/figure_3D_output.tab",quote = F, sep="\t")
write.list(figure_4Dl_output, file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/Figure3/figure_3E_output.tab",quote = F)
write.list(figure_4Dr_output, file = "/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/Figure3/figure_3F_output.tab",quote = F)


#Figure output
#All figure are saved in one file
pdf("/env/cns/home/nguerin/projet_CNM/Articles/PelagoNitro/Figure3/Figure_3.pdf", width = 9, height = 6)
print(fig3B_MAplot)
print(fig3C_MAplot)
print(fig3D_MAplot)
print(fig3E_euler)
print(fig3F_euler)
dev.off() 

```

