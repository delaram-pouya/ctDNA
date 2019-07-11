## in This script we'll try to run a pileup on the dataseries2 bam files, 
## then filter them based on vcf files and keep only the somatic positions,
## then remove the positions which are not present in all the files of each person(make them comparable)
## check correlation of VAFs in different files of each patient

## input : bam + vcf 
## output : VAF plots 

source('error/Codes/2_pileup_functions.R')
Initialize()

## loading the VCF files for P1 & P2

setwd('/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/Analysis/p1')
f1 <- list.files('.','lev2.txt$')
f1.vcf <- list.files('.','lev2.vcf$')
v1 <- mclapply(f1, fread, header=T, mc.cores = detectCores()-2)
v1.vcf <- mclapply(f1.vcf, read.vcfR, verbose = FALSE , mc.cores = detectCores()-2)
v1.vcf <- mclapply(v1.vcf, function(x)data.frame(getFIX(x)), mc.cores = detectCores()-2)



setwd('/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/Analysis/p2')
f2 <- list.files('.','lev2.txt$')
f2.vcf <- list.files('.','lev2.vcf$')
v2 <- mclapply(f2, fread, header=T, mc.cores = detectCores()-2)
v2.vcf <- mclapply(f2.vcf, read.vcfR, verbose = FALSE , mc.cores = detectCores()-2)
v2.vcf <- mclapply(v2.vcf, function(x)data.frame(getFIX(x)), mc.cores = detectCores()-2)



## cleaning and concatinating them in a list 
f <- gsub('_snv_lev2.txt','',c(f1,f2))
v <- c(v1,v2)
vcf <- c(v1.vcf, v2.vcf)
vcf <- lapply(vcf,function(x){x$CHROM<-as.character(x$CHROM);x$POS<-as.numeric(as.character(x$POS));
x$REF<-as.character(x$REF);x$ALT<-as.character(x$ALT);return(x[,!colnames(x)%in%c('ID','QUAL','FILTER')])})
names(v) <- f;names(vcf) <- f
total.vcf <- unique(do.call(rbind, vcf))



####### loading bam files
PATH = '/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/ALN/'

#### p1
l1 <- list.files(path=paste0(PATH,'p1'), pattern='*.addRG.bam$')
samples1 <- sub(".addRG.bam","",l1)
bf1 <- mclapply(l1, function(i) BamFile(i, index=paste0(i, ".bai")), mc.cores=detectCores()-2)

# bed.pileup1 <- mclapply(bf1, pileup, pileupParam=p_param, mc.cores=detectCores()-2)
bed.pileup1 <- readRDS('/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/ALN/p1/bedPileup_p1.rds')
names(bed.pileup1) <- samples1
sapply(1:length(bed.pileup1),function(i)bed.pileup1[[i]]$seqnames <<- as.character(bed.pileup1[[i]]$seqnames))

## filtreing the pileup on the position which are present in the vcf file 
### keeping only the somatic mutations

p1 <- sapply(1:length(bed.pileup1), function(i) {
  tmp <- subset(bed.pileup1[[i]],bed.pileup1[[i]]$seqnames %in% total.vcf$CHROM & 
                  bed.pileup1[[i]]$pos %in% total.vcf$POS);return(tmp)},simplify = F)

#saveRDS(p1,'/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/ALN/p1/p1.rds')



#### p2
l2 <- list.files(path=paste0(PATH,'p2'), pattern='*.addRG.bam$')
samples2 <- sub(".addRG.bam","",l2)
bf2 <- mclapply(l2, function(i) BamFile(i, index=paste0(i, ".bai")), mc.cores=detectCores()-2)
# bed.pileup2 <- mclapply(bf2, pileup, pileupParam=p_param, mc.cores=detectCores()-2)
bed.pileup2<- readRDS('/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/ALN/p2/bedpileup2_all.rds')
names(bed.pileup2) <- samples2

sapply(1:length(bed.pileup2),function(i)bed.pileup2[[i]]$seqnames <<- as.character(bed.pileup2[[i]]$seqnames))
p2 <- sapply(1:length(bed.pileup2), function(i) {
  tmp <- subset(bed.pileup2[[i]],bed.pileup2[[i]]$seqnames %in% total.vcf$CHROM & 
                  bed.pileup2[[i]]$pos %in% total.vcf$POS);return(tmp)},simplify = F)

#saveRDS(p2,'/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/ALN/p2/p2.rds')




##### generating p1 + p2 > ram issues
bed.pileup <- c(bed.pileup1, bed.pileup2)
samples <- c(samples1, samples2)
names(bed.pileup) <- samples

sapply(1:length(bed.pileup),function(i)bed.pileup[[i]]$seqnames <<- as.character(bed.pileup[[i]]$seqnames))
p <- sapply(1:length(bed.pileup), function(i) {
  tmp <- subset(bed.pileup[[i]],bed.pileup[[i]]$seqnames %in% vcf[[i]]$CHROM & bed.pileup[[i]]$pos %in% vcf[[i]]$POS)
  return(tmp)},simplify = F)

rm(bed.pileup1); rm(bed.pileup2)




##### loading p1 + p2 
p1 <- readRDS('error/Data/P1/p1.rds')
p2 <- readRDS('error/Data/P2/p2.rds')

## keeping only the position which are present in all files of a patient
p1 <- find.common(p1); names(p1) <- samples1
p2 <- find.common(p2); names(p2) <- samples2
p <- c(p1, p2)

freq <- mclapply(p, MakepileupFreq, mc.cores = detectCores()-2)
freq2 <- mclapply(freq, Add.attrib, mc.cores = detectCores()-2)
freq2 <- mclapply(freq2, function(x){x[is.na(x)]<-0; x$E<-(1-(x$F1+x$F2));return(x)})
names(freq2) <- samples

freq.tmp <- freq2[c(2,3,4,6,7,8)]
freq3 <- sapply(1:length(vcf),
                function(i) 
                  merge(freq.tmp[[i]],vcf[[i]],by.x=c('seqnames','start'),by.y=c('CHROM','POS'),all.x=T,all.y=F),
                simplify=F)

#freq3 <- lapply(freq3, subset, depth>0)
#freq3 <- lapply(freq3, subset, ALT %in% c('A','T','C','G') | is.na(ALT))


### calculating VAF
freq3 <- mclapply(freq3, function(x){
  x$VAF<-sapply(1:nrow(x), function(i){
    find.VAF(i,x)});x},mc.cores=detectCores()-2)

names(freq3) <- names(freq.tmp)
#saveRDS(freq3, 'simulator/Tumor_origin_read_probability/pileUpVAFall.rds')
freq3 <- readRDS('simulator/Tumor_origin_read_probability/pileUpVAFall.rds')
saveRDS(vcf, 'simulator/Tumor_origin_read_probability/allVCFs.rds')




#### Removing heterozygous positions
heter.p1 <- subset(freq2[[1]], F1<0.7 & F2>0.3); heter.p1 <- paste0(heter.p1$seqnames,sep='_',heter.p1$start)
heter.p2 <- subset(freq2[[5]], F1<0.7 & F2>0.3); heter.p2 <- paste0(heter.p2$seqnames,sep='_',heter.p2$start)

index <- sapply(1:length(freq3), function(i) {
  if(i %in% c(1,2,3)){paste0(freq3[[i]]$seqnames,sep='_',freq3[[i]]$start) %in% heter.p1}
  else if(i %in% c(4,5,6)){paste0(freq3[[i]]$seqnames,sep='_',freq3[[i]]$start) %in% heter.p2}
}) 

sapply(1:length(freq3), function(i){freq3[[i]]$homo <<- !index[[i]]},simplify = F)





########### All positions
#### p1
m1 <- sapply(1:3,function(i)
  merge(freq3[[1]], freq3[[i]], by.x=c('seqnames', 'start'),by.y=c('seqnames','start'),all.x=T,all.y=T),
  simplify=F)

m2 <- sapply(2:3,function(i)
  merge(freq3[[2]], freq3[[i]], by.x=c('seqnames', 'start'),by.y=c('seqnames','start'),all.x=T,all.y=T),
  simplify=F)


#### p2
m1 <- sapply(4:6,function(i)
  merge(freq3[[4]], freq3[[i]], by.x=c('seqnames', 'start'),by.y=c('seqnames','start'),all.x=T,all.y=T),
  simplify=F)

m2 <- sapply(5:6,function(i)
  merge(freq3[[5]], freq3[[i]], by.x=c('seqnames', 'start'),by.y=c('seqnames','start'),all.x=T,all.y=T),
  simplify=F)



m1.plot <- sapply(1:length(m1),function(i)ggplot(m1[[i]],aes(x=VAF.x,y=VAF.y,color=F1.y))+geom_point()+
         xlab('cfDNA1')+ylab(c('cfDNA1','cfDNA2','tissue')[i])+theme_minimal()+ggtitle('p2-VAF')+
         scale_color_gradientn(colours=rainbow(4)),simplify = F)

m2.plot <- sapply(1:length(m2),function(i)ggplot(m2[[i]],aes(x=VAF.x,y=VAF.y,color=F1.y))+geom_point()+
                    xlab('cfDNA2')+ylab(c('cfDNA2','tissue')[i])+theme_minimal()+ggtitle('p2-VAF')+
                    scale_color_gradientn(colours=rainbow(4)),simplify = F)


pdf('Error/Plots/VAF.pdf', width=4.7, height= 6)
m1.plot
m2.plot
dev.off()




######## Homo positions

m1.plot <- sapply(1:length(m1),function(i)
  ggplot(m1[[i]],aes(x=VAF.x,y=VAF.y,color=homo.x))+
    geom_point()+xlab('cfDNA1')+ylab(c('cfDNA1','cfDNA2','tissue')[i])+
    theme_minimal()+ggtitle('p2-VAF'),simplify = F)

m2.plot <- sapply(1:length(m2),function(i)
  ggplot(m2[[i]],aes(x=VAF.x,y=VAF.y,color=homo.x))+
    geom_point()+xlab('cfDNA2')+ylab(c('cfDNA2','tissue')[i])+
    theme_minimal()+ggtitle('p2-VAF'),simplify = F)


pdf('Error/Plots/VAF-homo.pdf', width=4.7, height= 6)
m1.plot
m2.plot
dev.off()




#############################

sapply(1:length(freq3), function(i)
  freq3[[i]]$sample<<-c('p1_cfDNA1','p1_cfDNA2','p1_tissue','p2_cfDNA1','p2_cfDNA2','p2_tissue')[i])

total <- do.call(rbind, freq3)

###############
sapply(1:length(freq2), function(i)freq2[[i]]$sample<<-samples[i])

Index <- sapply(1:length(freq2), function(i) {
  if(i %in% c(1,2,3,4)){paste0(freq2[[i]]$seqnames,sep='_',freq2[[i]]$start) %in% heter.p1}
  else if(i %in% c(5,6,7,8)){paste0(freq2[[i]]$seqnames,sep='_',freq2[[i]]$start) %in% heter.p2}
}) 

sapply(1:length(freq2), function(i){freq2[[i]]$homo <<- !Index[[i]]},simplify = F)
total2 <- do.call(rbind, freq2)





pdf('Error/Plots/error.pdf')

ggplot(subset(total,depth>20& homo& !is.na(ALT) ), aes(E+1e-7,color=sample))+
  geom_density()+scale_x_log10()+theme_minimal()+
  xlab('E+1e-07(log)')+ggtitle('error-distribution (lev2-depth>20-homo-SNVs)')


ggplot(subset(total,depth>20& homo& !is.na(ALT) ), aes(y=E+1e-7,x=sample))+
  geom_violin(width=0.7,aes(fill=sample))+scale_y_log10()+theme_minimal()+
  ylab('E+1e-07(log)')+ggtitle('error-distribution (lev2-depth>20-homo-SNVs)')


ggplot(subset(total2,depth>20& homo), aes(E+1e-7,color=sample))+
  geom_density()+scale_x_log10()+theme_minimal()+xlab('E+1e-07(log)')+
  ggtitle('error-distribution (lev2-depth>20-homo-All)')

ggplot(subset(total2,depth>20& homo), aes(y=E+1e-7,x=sample))+
  geom_violin(width=0.7,aes(fill=sample))+scale_y_log10()+theme_minimal()+
  ylab('E+1e-07(log)')+ggtitle('error-distribution (lev2-depth>20-homo-All)')


dev.off()




pdf('Error/Plots/VAF-dis.pdf')

ggplot(subset(total, depth>20& homo & !is.na(ALT)), aes(VAF+1e-7,color=sample))+
  geom_density()+scale_x_log10()+theme_minimal()+ggtitle('VAF-distribution (lev2-d>20-homo-SNVs)')

ggplot(subset(total, depth>20& homo), aes(VAF+1e-7,color=sample))+
  geom_density()+scale_x_log10()+theme_minimal()+ggtitle('(lev2-d>20-homo-All)')

ggplot(subset(total, depth>20& homo & !is.na(ALT)), aes(VAF+1e-7,x=sample))+
  geom_violin(width=0.7,aes(fill=sample))+theme_minimal()+
  ggtitle('VAF-distribution (lev2-d>20-homo-SNVs)')+ylab('VAF.log')+scale_y_log10()

ggplot(subset(total, depth>20& homo ), aes(VAF+1e-7,x=sample))+
  geom_violin(width=0.7,aes(fill=sample))+theme_minimal()+
  ggtitle('VAF-distribution (lev2-d>20-homo-All)')+ylab('VAF.log')+scale_y_log10()

dev.off()



ggplot(subset(freq3[[1]],depth>20& homo & !is.na(ALT)),aes(x=F2,y=VAF))+geom_point()

