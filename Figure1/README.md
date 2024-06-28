# Figure 1
This directory contains the files and scripts to generate Figure 1 of the article.   
Input files are :  
Gene expression levels of P. calceolata RCC100 in Tara samples : https://zenodo.org/records/6983365  
Abundance of Pelagomonas calceolata RCC100 in Tara Oceans samples : PelagoV4_metaG_relative_abundance.tab  in this directory
Metadata of Tara Oceans expedition : Environmental_metadata.tab in this directory  
Coordinates of Tara Oceans stations : TaraCoord.txt in this direcory  
Coordinates of Tara Pacific stations : Coordinates_islands_TO-TP.txt 

The following code was executed on R version 4.1.1

## Library loading
```{r}
library(ggplot2)
library(tidyr) 
library(dplyr)
library(ggforce)
library(viridis)
library(ggrepel)
library(DESeq2)
library(tibble)
```

## Figure 1A
Nitrate concentrations measured during _Tara_ Oceans expedition. The colour code indicates nitrate concentrations in µmol/l for surface and DCM samples in the upper and bottom part of each dot respectively.

### Inputs and processing
```{r}
#loading Tara Oceans metadata
metadata <- read.table("Environmental_metadata.tab", row.names=1)
tab <-  subset(metadata, select = -c(Latitude, Longitude))
#Extracting station Names
row.names(tab) -> tab$Names
#row.names(tab) -> tab$Sample_bis

#loading map data
mapWorld <- map_data('world')

#loading coordinates and arrows from Tara Oceans
coord<-read.table("TaraCoord.txt",sep="\t",h=T,stringsAsFactors = F)
#merging metadata and coordinates
tab <- merge(tab, coord[,c(1:5)], by.x="Sample_label", by.y="Stations", all = T)

#loading additionnal coordinates from Tara Oceans
coord<-read.table("Coordinates_islands_TO-TP.txt",sep="\t",h=T,stringsAsFactors = F)
#Extracting station numbers
tab$Station_number<- sub("\\D.*", "", tab$Sample_label)
#merging metadata and coordinates
tab <- merge(tab, coord[,c(1,3,4)], by.x="Station_number", by.y="ID")

#if xend is emppty then it takes Longitude value
ifelse(is.na(tab$xend)==T, tab$xend<-tab$Longitude, NA)
ifelse(is.na(tab$yend)==T, tab$yend<-tab$Latitude, NA)

#remove NA
tab <- tab %>% drop_na(Nitrate_median)

#if there is two size fraction, keep only the GGMM, else keep the GGZZ in order to have a map without duplicates
tab <- mutate(tab,
               Fraction = case_when(
                 grepl("GGMM", Names) ~ "GGMM",
                 grepl("GGZZ", Names) ~ "GGZZ"))
length(grep("GGMM", tab$Fraction)) #68
length(grep("GGZZ", tab$Fraction)) #44

#remove fraction duplicates
#if isUnique== F, then delete the GGZZ, else do nothing
tab <- tab %>%
  arrange(Sample_label, Fraction) %>%
  filter(duplicated(Sample_label) == FALSE)
length(grep("GGMM", tab$Fraction)) #67
length(grep("GGZZ", tab$Fraction)) #6

#creation of map_table, a copy of tab that will be used to generate the map (figure 1A)
figure_1A_output <- tab
#addition of half circle according to depth
figure_1A_output <- mutate(figure_1A_output,
       start = case_when(
         Layer=="DCM" ~ -pi/2,
         Layer=="SUR" ~ pi/2))
```

### Plot and Output
```{r}
## Figure 1A 
fig1A_mapNO3 <- ggplot()+
    geom_polygon(data=mapWorld, aes(x=long, y=lat, group = group),fill="grey70",color="grey70") + 
    geom_segment(data=figure_1A_output, aes(x=x,xend=xend,y=y,yend=yend),linewidth=0.1,arrow=arrow(length = unit(1, "mm")))+ 
    geom_arc_bar(data=figure_1A_output[is.na(figure_1A_output$x)==F,], aes(x0 = x, y0 = y, r0 = 0, r=4, start = start, end = start + pi, fill = Nitrate_median),color = NA) + # r = normalised_metaT
    coord_quickmap(xlim =c(-155,75), ylim = c(-70,50))+
    geom_arc_bar(data=figure_1A_output[is.na(figure_1A_output$x)==T,], aes(x0 = Longitude, y0 = Latitude, r0 = 0, r = 4, start = start, end = start + pi, fill = Nitrate_median),color = NA) +
    scale_fill_viridis(option= "C", trans="log10", na.value="#0D1687") +
    theme(legend.position=c(0.15,0.18),legend.justification=c(1,0.5),legend.background=element_rect(fill="white",color="black",size=0.1)) + theme_bw() +
    labs(x="", y="")
#figure 1A
fig1A_mapNO3
#output for figure 1A
figure_1A_output
```

## Figure 1B
Relative abundance of _P. calceolata_ in _Tara_ samples estimated from metagenomics reads according to the concentration of nitrate (µM).

### Inputs and processing
```{r}
#tab, generated previously
tab
#for tab, create Sample_bis to match Sample_bis column from relative_abundance
tab$Sample_bis <-gsub('.{2}$', '', tab$Names)

#loading metagenomic relative abundance data (from Guérin et al. 2020)
relative_abundance <- read.table("PelagoV4_metaG_relative_abundance.tab", header=T)

#for Tara Oceans samples, extracting sample label 
relative_abundance$Sample_label <- gsub('.{2}$', '', relative_abundance$Sample_bis)
#merge tab and relative_abundance into figure_1B_output
figure_1B_output <- merge(tab, relative_abundance, by="Sample_bis",all=F)

```
### Plot and Output
```{r}
fig1B_distri <- ggplot(data=figure_1B_output, aes(x=Nitrate_median, y=relative_abundance,color=Layer)) +
  geom_point()+
  scale_color_manual(values = c("#A4262C", "#1d3084ff"))+
  #scale_x_continuous(trans = 'log2')+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 2))+
  geom_vline(xintercept = 2, color="red")+
  theme_bw()+
  labs(x="Nitrate concentration (µM)", y="Metagenomic relative abundance")
#figure 1B
fig1B_distri
#output for figure 1B
figure_1B_output
```
## Figure 1C
_P. calceolata_ gene-expression levels between low-nitrate (NO3 < 2 µM, n=69) and high-nitrate samples (NO3 > 2 µM, n=43). Log2FC between low- and high-nitrate samples are given according to their mean expression level (normalized with DESeq2). Differentially expressed genes with p-value < 0.01 and log2FC >1 or < −1 are coloured in blue.

### Input and processing
```{r}
#loading raw metatranscriptomic count table
rawCountTable <- read.table("20230427_RCC100-Nitrate_transcriptomes_rawcounts.tsv", check.names=FALSE)

#loading metadata
sampleInfo <- metadata

#transform metadata to set the nitrate concentration threshold to 2 µM
sampleInfo <- mutate(sampleInfo, nitrate_quantity = case_when(
  Nitrate_median < 2 ~ "low",
  Nitrate_median > 2 ~ "high"))

#removes samples for which nitrate concentration data is not available 
row_names_to_remove<-rownames(sampleInfo[is.na(sampleInfo$nitrate_quantity)==T,])
sampleInfo <- sampleInfo[!(row.names(sampleInfo) %in% row_names_to_remove),]
rm(row_names_to_remove)

#remove samples from rawCountTable that are not found in sampleInfo
rawCountTable <- rawCountTable[,(colnames(rawCountTable) %in% row.names(sampleInfo))]
#ordering sampleInfo
sampleInfo <- sampleInfo[ order(row.names(sampleInfo)), ]

#We are left with 112 samples in each tables
dim(rawCountTable)
dim(sampleInfo)
dim(filter(sampleInfo, nitrate_quantity == "low")) #69 low NO3
dim(filter(sampleInfo, nitrate_quantity == "high")) #43 high NO3

#Creating DESeq2 matrix
se_star_matrix <- DESeqDataSetFromMatrix(countData = rawCountTable,
                                         colData = sampleInfo,
                                         design = ~ nitrate_quantity)
#Creating DESeq2 object
dds <- DESeq(se_star_matrix)

#log2 fold change calculation
res <- lfcShrink(dds, coef=resultsNames(dds)[2], type = 'normal')

#mantel plot can be generated
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, xlim=xlim, ylim=ylim, main="DEG according to environmental nitrate concentration", alpha=0.01)
abline(h=2, col="red")
abline(h=-2, col="red")

#table res contains the normalised results of the DEG analysis
#formatting res into figure_1C_output in order to create the figure 1C
figure_1C_output <- as.data.frame(res)
figure_1C_output <- mutate(figure_1C_output,
       regulation = case_when(
         log2FoldChange > 1 & padj<0.01 ~ "upregulated",
         log2FoldChange < -1 & padj<0.01 ~ "downregulated"))
#creation of figure_1C_output
top20_genes <- figure_1C_output %>%
  rownames_to_column('Gene') %>%
  filter(regulation=="upregulated" | regulation=="downregulated") %>%
  group_by(regulation) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  dplyr::slice(1:10)
```

### Plot and Output
```{r}
fig1C_MAplot <- ggplot()+  
  geom_point(data=figure_1C_output, aes(x=log2(baseMean), y=log2FoldChange, color=regulation), size=1)+  
  geom_text_repel(data=top20_genes, aes(x=log2(baseMean), y=log2FoldChange, label=Gene), size=3, point.size = NA) +
  ylim(-3,3) +
  geom_hline(yintercept = 0) + geom_hline(yintercept = c(-2, 2), color="red") +
  labs(title = "DEG according to environmental nitrate concentration", x="Mean of normalized count (TPM)", y="log2FC") +
  theme_bw() + 
  scale_color_manual(values = c("#0000b5", "#00bef2"), na.value = "grey80",name = "padj"
                     , labels = c("overexpressed in low nitrate","underexpressed in low nitrate",  "non significative"))

fig1C_MAplot
figure_1C_output
```
