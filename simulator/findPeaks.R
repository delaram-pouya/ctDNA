## In this script, we'll check if the are any regions that 
## lead to high depth in all cfDNA files(to be used in targeted cfDNA panels)  
## bedtools converage has been applied to the data, it tries 
## to find the mean depth for each region of a given bed file

Initialize <- function(){
  library(plyr)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(ggplot2)
  library(parallel)
}

Initialize()
setwd('/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/Raw_files/')


f <- list.files('.','.bedgraph$')
f <- f[-c(1,2,9,10)]
d <- mclapply(f, fread, header=F, mc.cores = detectCores()-2)
names <- sub('MeanCovBed_','', sub('.bedgraph','',f))
names <- sub('_sort','',names)
d <- lapply(d, function(x){colnames(x)=c('chr','start','end','mean');x})
d <- lapply(d, function(x){x$width=x$end-x$start;x})
names(d) <- names
depth <- lapply(d,function(x) x$mean)
names(depth) <- names
df<- data.frame(do.call(cbind,depth))

colnames(df) = c('blood','blood','cfDNA','cfDNA','cfDNA','cfDNA','tissue','tissue')
df_melted <- melt(df)
colnames(df_melted) = c('sample', 'value')

df_melted$sample <- as.character(df_melted$sample)
ggplot(df_melted, aes(x=value))+geom_density(aes(fill=sample),alpha=0.3)+
  scale_x_log10()+theme_bw()+ggtitle('mean depth distribution')+xlab('mean depth(log-transformed)')


##>>> all high depth regions in tissue lead to high depth regions in cfDNA
regions  <- do.call(cbind, d)
cfDNAs = subset(regions, cfDNA1_p1.mean<30 & cfDNA1_p2.mean<30 & cfDNA2_p1.mean<30 & 
                  cfDNA2_p2.mean<30 & tissue_p1.mean>30 & tissue_p2.mean>30 & blood_p1.mean>30 & blood_p2.mean>30)
cfDNAs = cfDNAs[,c('cfDNA1_p1.chr','cfDNA1_p1.start','cfDNA1_p1.end','cfDNA1_p1.mean','cfDNA1_p2.mean',
                   'cfDNA2_p1.mean','cfDNA2_p2.mean','tissue_p1.mean','tissue_p2.mean', 'blood_p1.mean','blood_p2.mean')]
cfDNAs


#### there are some high-depth regions in cfDNA which have low depth in tissue >>> all in chr11
tis = subset(regions, cfDNA1_p1.mean>100 & cfDNA1_p2.mean>100 & cfDNA2_p1.mean>100 & 
                  cfDNA2_p2.mean>100 & tissue_p1.mean<30 & tissue_p2.mean<30 & blood_p1.mean<30 & blood_p2.mean<30)
tis = tis[,c('cfDNA1_p1.chr','cfDNA1_p1.start','cfDNA1_p1.end','cfDNA1_p1.mean','cfDNA1_p2.mean',
                   'cfDNA2_p1.mean','cfDNA2_p2.mean','tissue_p1.mean','tissue_p2.mean', 'blood_p1.mean','blood_p2.mean')]
tis






##### Visualization ####

pdf('depthPlots.pdf',width=13)
#### P1
p1=ggplot(df,aes(x=blood_p1,y=cfDNA1_p1))+geom_point()+theme_bw()+ggtitle('P1')
p2=ggplot(df,aes(x=blood_p1,y=cfDNA2_p1))+geom_point()+theme_bw()
p3=ggplot(df,aes(x=tissue_p1,y=cfDNA1_p1))+geom_point()+theme_bw()
p4=ggplot(df,aes(x=tissue_p1,y=cfDNA2_p1))+geom_point()+theme_bw()
p5=ggplot(df,aes(x=cfDNA1_p1,y=cfDNA2_p1))+geom_point()+theme_bw()
p6=ggplot(df,aes(x=blood_p1,y=tissue_p1))+geom_point()+theme_bw()
gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2,ncol=3)

#### P2
p1=ggplot(df,aes(x=blood_p2,y=cfDNA1_p2))+geom_point()+theme_bw()+ggtitle('P1')
p2=ggplot(df,aes(x=blood_p2,y=cfDNA2_p2))+geom_point()+theme_bw()
p3=ggplot(df,aes(x=tissue_p2,y=cfDNA1_p2))+geom_point()+theme_bw()
p4=ggplot(df,aes(x=tissue_p2,y=cfDNA2_p2))+geom_point()+theme_bw()
p5=ggplot(df,aes(x=cfDNA1_p2,y=cfDNA2_p2))+geom_point()+theme_bw()
p6=ggplot(df,aes(x=blood_p2,y=tissue_p2))+geom_point()+theme_bw()
gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2,ncol=3)

##### P1 vs P2
p1=ggplot(df,aes(x=blood_p1,y=blood_p2))+geom_point()+theme_bw()+ggtitle('P1 vs P2')
p2=ggplot(df,aes(x=tissue_p1,y=tissue_p2))+geom_point()+theme_bw()
p3=ggplot(df,aes(x=cfDNA1_p1,y=cfDNA1_p2))+geom_point()+theme_bw()
p4=ggplot(df,aes(x=cfDNA2_p1,y=cfDNA2_p2))+geom_point()+theme_bw()
gridExtra::grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)

