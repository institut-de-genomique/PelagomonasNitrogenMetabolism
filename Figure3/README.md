# Figure 3  
This directory contains the files and scripts to generate Figure 3 of the article.  
Input files are :  
Cell count for P. calceolata _RCC100_ in 3 different conditions in triplicates : GrowthCyanateAxenic.tsv  in this directory

The following code was executed on R version 4.1.1

## Required libraries  
```r
library(reshape2)
library(ggplot2)
library(growthcurver)
```
## Data loading and processing
```r
tab<-read.table("GrowthCyanateAxenic.tsv",sep="\t",h=T,stringsAsFactors = FALSE)
#Transformation long format
data<-melt(tab,measure.vars = c(2:ncol(tab)),variable.name = "Replicates")
#Ajout d'une colonne pour conditions (hors réplicats)
data$Condition<-sub("R.*","",data$Replicates)
#Min,Max,Mean per Condition
data2<-do.call(data.frame,aggregate(data,value~Condition+time,function(x){c(Min=min(x),Max=max(x),Mean=mean(x))}))
#format require by GrowthCurve package (short format and column time)
data3<-data.frame(acast(data2,time~Condition,value.var = "value.Mean"))
data3$time<-as.numeric(rownames(data3))
```
## Fit with GrowthCurver package
```r
Fit <-SummarizeGrowthByPlate(data3,bg_correct = "none")
write.table(Fit,file="FitMetrics_Average.tab",quote=F,sep="\t",row.names = F)
#Curve values based on Fit metrics
FitCurve<-do.call(rbind.data.frame,apply(Fit,1,function(x){
  data.frame(sample=rep(x[1],max(data3$time)+5),
             Count=as.numeric(x[2])/(1+((as.numeric(x[2])-as.numeric(x[3]))/as.numeric(x[3]))*exp(-as.numeric(x[4])*c(1:(max(data3$time)+5)))),
             time=c(1:(max(data3$time)+5)))
  }))
```
## Plot
```r
pdf(file="GrowthCurveFit.pdf")
ggplot(FitCurve)+
  geom_line(aes(x=time,y=Count,color=sample),linewidth=1.5)+
  geom_point(data=Fit[Fit$r>0.1,],aes(x=t_mid,y=k/2,color=sample),shape=18,size=5)+
  geom_point(data=data2,aes(x=time,y=value.Mean),color="black")+
  geom_errorbar(data=data2,aes(ymax=value.Max,ymin=value.Min,x=time),size=1)+
  geom_label(data=Fit,aes(x=0,
                             y=seq(max(data$value)/2,max(data$value),length.out=nrow(Fit)),
                             label=paste(sample,": r=",round(r,digits = 2),sep=""),
                             color=sample),
             hjust=0)+
  theme(legend.position="none",panel.background = element_blank(),panel.grid=element_line(colour="grey75",linetype=2))+
  scale_color_brewer(palette ="Set1")
```
