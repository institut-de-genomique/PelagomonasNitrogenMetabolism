#Supplementary Figure 1
This directory contains the files and scripts to generate Supplementary Figure 1 (panels C and D) of the article.
Input files are :  
Chlorophyll a fluorescence for RCC100 and RCC697 in different culture conditions:  Data-SupFig1.tab  

The following code is executed on R version 4.1.1  
```{r}

tab<-read.table("Data-SupFig1.tab",h=T,sep="\t")
#Miniman, maximum and average for each condition, strain and time point
tab2<-do.call(data.frame,aggregate(tab,Fluorescence~Condition+Day+Strain,function(x){c(Min=min(x),Max=max(x),Mean=mean(x))}))
tab2$Condition<-paste(tab2$Condition,"µM")
#plot for panels C and D
pdf(file="GrowthCurvesLowN.pdf",height=9)
ggplot(tab2)+
  geom_line(aes(x=Day,y=Fluorescence.Mean,color=Condition))+
  geom_errorbar(aes(ymax=Fluorescence.Max,ymin=Fluorescence.Min,x=Day),size=0.5,width=0.2)+
  geom_point(aes(x=Day,y=Fluorescence.Mean,color=Condition))+
  theme(legend.position="bottom",panel.background = element_blank(),panel.grid=element_line(colour="grey75",linetype=2))+
  scale_color_brewer(palette ="Set1")+
  scale_y_continuous(trans="log10")+
  facet_grid(Strain~.,scales = "free")
dev.off()
```
