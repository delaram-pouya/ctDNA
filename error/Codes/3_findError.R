## In this script, we'll try to find the general and base-wise error distribution in each sample 
## first, we'll find the homozygous positions based on the blood samples(positions with F1>0.8) 
## the, we'll remove the heterozygous positions from all samples and 
## remove somatic mutations(based on the VCF files) from the tissue and cfDNA samples 
## calculate error rate for each position as E= F2+F3+F4 (hetero and somatic positions have been removed)
## make the pair-wise error table

## input: pileup files 
## output: base-wise error distribution plots for each sample seperately

source('Error/Codes/4_findError_functions.R')
Initialize()


#### Importing all pileup files

bed.pileup1 <- readRDS('Error/Data/P1//bedPileup_p1.rds')
bed.pileup2 <- readRDS('Error/Data/P2/bedpileup2_all.rds')
samples <- c('p1_blood','p1_cfDNA1','p1_cfDNA2','p1_tissue',
             'p2_blood','p2_cfDNA1','p2_cfDNA2','p2_tissue')


## HOMOZYGOUS 
#### finding homozygous positions for p1 & p2 based on blood sample

blood <- list(bed.pileup1[[1]],bed.pileup2[[1]])
blood <- mclapply(blood, data.table, mc.cores = detectCores()-2)
bloodPileup <- mclapply(blood, preprop, mc.cores = detectCores()-2)
bloodPileup <- mclapply(bloodPileup, MakepileupFreq.eff, mc.cores = detectCores()-2)
bloodPileup <- mclapply(bloodPileup, Add.depth , mc.cores = detectCores()-2)
bloodPileup2 <- mclapply(bloodPileup, Add.max , mc.cores = detectCores()-2) # to remove hetero: Add.max.blood()
names(bloodPileup2) <- c('p1_blood','p2_blood')
homo <- mclapply(bloodPileup2, function(x){df=subset(x,F1>0.8);paste0(x$seqnames,sep='_',x$pos)})
homo <- readRDS('Error/Data/homozygous_pos.rds')





######### TISSUE ##############
#### Somatic mutations

v1 <- fread('Error/Data/P1/VCFs/p1_tissue_snv_lev3.txt', header = T)
v2 <- fread('Error/Data/P2/VCFs/p2_tissue_snv_lev3.txt', header = T)
v <- list(v1,v2); names(v)=c('p1_tissue','p2_tissue')
somatics <- mclapply(v, function(x)paste0(x$chrom,sep='_',x$pos),mc.cores = detectCores()-2)


#### importing tissue pileup files 
tissue <- list(bed.pileup1[[4]],bed.pileup2[[4]])
tissue <-  mclapply(tissue, data.table, mc.cores = detectCores()-2)



# remove hetero positions
index_homo <- sapply(1:length(tissue), 
                     function(i) paste0(tissue[[i]]$seqnames,sep='_',tissue[[i]]$pos) %in% homo[[i]],
                     simplify=F)

tissueHomo <- sapply( 1:length(tissue), function(i) tissue[[i]][ index_homo[[i]] ,],simplify = F )
saveRDS(tissueHomo, 'Error/Data/tisPileupHomoFilt.rds')



# remove somatic mutations
index_somatic <- sapply(1:length(tissueHomo), 
                        function(i) paste0(tissueHomo[[i]]$seqnames,sep='_',tissueHomo[[i]]$pos) %in% somatics[[i]],
                        simplify=F)

tissueHomoNoMut <- sapply( 1:length(tissueHomo), function(i)tissueHomo[[i]][ (!index_somatic[[i]]) ,],simplify = F )

tissueHomoNoMut <- mclapply(tissueHomoNoMut, preprop, mc.cores = detectCores()-2)
tissueHomoNoMut <- mclapply(tissueHomoNoMut, MakepileupFreq.eff, mc.cores = detectCores()-2)
tissueHomoNoMut <- mclapply(tissueHomoNoMut, Add.depth , mc.cores = detectCores()-2)
tissueHomoNoMut <- mclapply(tissueHomoNoMut, Add.max , mc.cores = detectCores()-2) 
names(tissueHomoNoMut) <- c('p1_tissue','p2_tissue') 
saveRDS(tissueHomoNoMut, 'Error/Data/tisPileupFinal.rds')



####  TISSUE error distributions ###
## Base-wise ##

# making the error table
bases <- c('A','T','C','G')
mat <- data.table(matrix(list(numeric(1)),ncol=4,nrow=4))
colnames(mat) <- bases; rownames(mat) <- bases
#tissueMat1 <- matUpdate(mat,subset(a2[[1]], select=c('A','T','C','G','m1','depth')))
#tissueMat2 <- matUpdate(mat,subset(a2[[2]], select=c('A','T','C','G','m1','depth')))
tissueMat1 <- readRDS('Error/Data/ErrorMatrix/tissueMat1.rds')
tissueMat2 <- readRDS('Error/Data/ErrorMatrix/tissueMat2.rds')

## visualization
makePlot(makeErrorTable(tissueMat1), 'tissue_p1')
makePlot(makeErrorTable(tissueMat2), 'tissue_p2')





####  BLOOD error distributions ###

## Run Add.max.blood function on pileup to remove the heterozygous positions
bloodPileupHomo <- mclapply(bloodPileup, Add.max.blood , mc.cores = detectCores()-2)
#saveRDS(bloodPileupHomo,'ALN/p1/bloodfinalpileup.rds')

### Total Error ###
## homo-blood:  F2+F3+F4 = 1-F1 = PCR error + Seq error (baseliene)

bloodPileupHomo[[1]]$sample <- 'p1_blood'; bloodPileupHomo[[2]]$sample <- 'p2_blood'
bloodPileupHomo.melt <- do.call(rbind, bloodPileupHomo)

p1=ggplot(bloodPileupHomo.melt,aes(x=ALT.freq,color=sample))+geom_density()+
  theme_bw()+xlab('F2+F3+F4')+ggtitle('seq error')

p2=ggplot(bloodPileupHomo.melt,aes(x=ALT.freq+1e-9,color=sample))+geom_density()+
  scale_x_log10()+theme_bw()+xlab('log(F2+F3+F4 + 1e-9)')+ggtitle('log seq error')

p3=ggplot(subset(bloodPileupHomo.melt, ALT.freq>0),aes(x=ALT.freq,color=sample))+
  geom_density()+theme_bw()+xlab('F2+F3+F4')+ggtitle('Non-zero seq error')

p4=ggplot(subset(bloodPileupHomo.melt, ALT.freq>0),aes(x=ALT.freq,color=sample))+scale_x_log10()+
  geom_density()+theme_bw()+xlab('log(F2+F3+F4)')+ggtitle('log Non-zero seq error')

zeroInfl = data.frame(table(bloodPileupHomo.melt$ALT.freq==0)); colnames(zeroInfl)=c('variable','Freq')
zeroInfl$variable = c('non-Zero','Zero') 
zeroInfl$Freq = round(zeroInfl$Freq/sum(zeroInfl$Freq),2)
p5=ggplot(zeroInfl,aes(x=variable,y=Freq,fill=variable))+ geom_bar(stat="identity",color='black',width=0.5)+theme_bw()

pdf('error/Plots/BloodErrorTotalMultPlots.pdf',height=12,width=10)
grid.arrange(p1,p2,p3,p4,p5,nrow=3,ncol=2)
descdist(bloodPileupHomo.melt$ALT.freq)
descdist(subset(bloodPileupHomo.melt, ALT.freq>0)$ALT.freq)
dev.off()




### Base-wise ###

# make the error table
bases <- c('A','T','C','G')
mat <- data.table(matrix(list(numeric(1)),ncol=4,nrow=4))
colnames(mat) <- bases; rownames(mat) <- bases
bloodMat1 <- matUpdate(mat,subset(bloodPileupHomo[[1]], select=c('A','T','C','G','m1','depth')))
bloodMat2 <- matUpdate(mat,subset(bloodPileupHomo[[2]], select=c('A','T','C','G','m1','depth')))
bloodMat1 <- readRDS('ALN/depth_all/ErrorMatrix/bloodMat1.rds')
bloodMat2 <- readRDS('ALN/depth_all/ErrorMatrix/bloodMat2.rds')

## visualization
makePlot(makeErrorTable(bloodMat1),'blood_p1')
makePlot(makeErrorTable(bloodMat2),'blood_p2')









####  cfDNA error distributions ###

#### finding Somatic mutations

cfdna1VCF <- lapply(list.files('Error/Data/P1/VCFs/','lev3.txt$', full.names = T)[1:2], read.table, header=T)
cfdna2VCF <- lapply(list.files('Error/Data/P2/VCFs/','lev3.txt$', full.names = T)[1:2], read.table, header=T)
cfDNAvcf = c(cfdna1VCF,cfdna2VCF)
names(cfDNAvcf)=c('p1_cfDNA1','p1_cfDNA2','p2_cfDNA1','p2_cfDNA2')
cfDNAvcf <- lapply(cfDNAvcf, function(x){colnames(x)=c('chrom','pos'); return(x[2:nrow(x),]) })
somatics.cfDNA <- mclapply(cfDNAvcf, function(x)paste0(x$chrom,sep='_',x$pos),mc.cores = detectCores()-2)


#### loading cfDNA pileups 

cfdnaPileup1 <- list(bed.pileup1[[2]], bed.pileup1[[3]])
cfdnsPileup2 <- list(bed.pileup2[[2]], bed.pileup2[[3]])
cfdnaPileup <- list(cfdnaPileup1[[1]], cfdnaPileup1[[2]], cfdnaPileup2[[1]], cfdnaPileup2[[2]])
cfdnaPileup.table <-  lapply(cfdnaPileup, data.table)


# remove hetero positions

cfHomo <- list(homo[[1]],homo[[1]],homo[[2]],homo[[2]]) # first two for P1, last two for P2
index_homo <- sapply(1:length(cfdnaPileup), function(i) 
  paste0(cfdnaPileup[[i]]$seqnames,sep='_',cfdnaPileup[[i]]$pos) %in% cfHomo[[i]],simplify=F)

cfdnaHomo <- sapply( 1:length(cfdnaPileup), function(i) cfdnaPileup[[i]][ index_homo[[i]] ,],simplify = F)
saveRDS(cfdnaHomo, 'Error/Data/cfdnaPileupHomoFilt.rds')
cfdnaHomo <- readRDS('Error/Data/cfdnaPileupHomoFilt.rds')


# remove somatic mutations
index_somatic <- sapply(1:length(cfdnaHomo), 
                        function(i) paste0(cfdnaHomo[[i]]$seqnames,sep='_',cfdnaHomo[[i]]$pos) %in% somatics.cfDNA[[i]],
                        simplify=F)

cfdnaHomoNoMut <- sapply(1:length(cfdnaHomo), function(i) cfdnaHomo[[i]][ (!index_somatic[[i]]) ,],simplify = F )

cfdnaHomoNoMut <- mclapply(cfdnaHomoNoMut, preprop, mc.cores = detectCores()-2)
cfdnaHomoNoMut <- mclapply(cfdnaHomoNoMut, MakepileupFreq.eff, mc.cores = detectCores()-2)
cfdnaHomoNoMut <- mclapply(cfdnaHomoNoMut, Add.depth , mc.cores = detectCores()-2)
cfdnaHomoNoMut <- mclapply(cfdnaHomoNoMut, Add.max , mc.cores = detectCores()-2) 
names(cfdnaHomoNoMut) <- c('p1_cfDNA1','p1_cfDNA2','p2_cfDNA1','p2_cfDNA2')
saveRDS(cfdnaHomoNoMut, '/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/Raw_files/cfDNAsPileupFinal.rds')


####  cfDNA error distributions ###
## Base-wise ##

# making the error table
mat <- data.table(matrix(list(numeric(1)),ncol=4,nrow=4))
colnames(mat) <- bases; rownames(mat) <- bases
#cfDNAMat1 <- matUpdate(mat,subset(cfdnaHomoNoMut[[1]], select=c('A','T','C','G','m1','depth')))
cfDNAMat1 <- readRDS('Error/Data/ErrorMatrix/cfDNAMat1.rds')
cfDNAMat2 <- readRDS('Error/Data/ErrorMatrix/cfDNAMat2.rds')
cfDNAMat3 <- readRDS('Error/Data/ErrorMatrix/cfDNAMat3.rds')
cfDNAMat4 <- readRDS('Error/Data/ErrorMatrix/cfDNAMat4.rds')
names(cfdnaHomoNoMut) <- c('p1_cfDNA1','p1_cfDNA2','p2_cfDNA1','p2_cfDNA2')

## visualization
makePlot(makeErrorTable(cfDNAMat1), 'p1_cfDNA1')
makePlot(makeErrorTable(cfDNAMat2), 'p1_cfDNA2')
makePlot(makeErrorTable(cfDNAMat3), 'p2_cfDNA1')
makePlot(makeErrorTable(cfDNAMat4), 'p2_cfDNA2')






######################################################## making a list of error dataframes
baseWiseErrorBlood1 <- GetBaseWiseError(bloodMat1)
baseWiseErrorBlood2 <- GetBaseWiseError(bloodMat2)

baseWiseErrorTissue1 <- GetBaseWiseError(tissueMat1)
baseWiseErrorTissue2 <- GetBaseWiseError(tissueMat2)

baseWiseError_cfDNA1 <- GetBaseWiseError(cfDNAMat1)
baseWiseError_cfDNA2 <- GetBaseWiseError(cfDNAMat2)
baseWiseError_cfDNA3 <- GetBaseWiseError(cfDNAMat3)
baseWiseError_cfDNA4 <- GetBaseWiseError(cfDNAMat4)



ListOfBaseErrorAllSamples = sapply(1:length(baseWiseErrorBlood2),
                                   function(i) {rbind(
                                     data.frame(val=baseWiseErrorTissue1[[i]],sample='tissue1'),
                                     data.frame(val=baseWiseErrorTissue2[[i]],sample='tissue2'),
                                     
                                     data.frame(val=baseWiseErrorBlood1[[i]],sample='blood1'),
                                     data.frame(val=baseWiseErrorBlood2[[i]],sample='blood2'), 
                                     
                                     data.frame(val=baseWiseError_cfDNA1[[i]],sample='p1_cfDNA1'), 
                                     data.frame(val=baseWiseError_cfDNA2[[i]],sample='p1_cfDNA2'), 
                                     data.frame(val=baseWiseError_cfDNA3[[i]],sample='p2_cfDNA1'), 
                                     data.frame(val=baseWiseError_cfDNA4[[i]],sample='p2_cfDNA2') 
                                     )}
                                   ,simplify = F)

names(ListOfBaseErrorAllSamples) = names(baseWiseErrorBlood1)



##################### compare different samples for each error  

#### including Zeros

pdf('Analysis/errorAnalysis/errorAnalysisPlots/compareSamples2/ListOfBaseErrorAllSamplesZeroIncluded2.pdf',width=14,height=7)
sapply(1:length(ListOfBaseErrorAllSamples),function(i){
  
  p1=ggplot(ListOfBaseErrorAllSamples[[i]],aes(val+1e-7,color=sample))+
    geom_density()+scale_x_log10()+ggtitle(names(ListOfBaseErrorAllSamples)[i])+theme_bw();
  
  p2=ggplot(ListOfBaseErrorAllSamples[[i]],aes(y=val+1e-7,x=sample))+
    geom_violin(aes(fill=sample))+scale_y_log10()+ggtitle(names(ListOfBaseErrorAllSamples)[i])+theme_bw();
  
  gridExtra::grid.arrange(p1,p2,nrow=1,ncol=2)
  },simplify = F )
dev.off()




######## Non-Zero errors

pdf('Analysis/errorAnalysis/errorAnalysisPlots/compareSamples2/ListOfBaseErrorAllSamplesNonZero2.pdf',width=23,height = 6)
sapply(1:length(ListOfBaseErrorAllSamples),function(i){
  
  BaseError = subset(ListOfBaseErrorAllSamples[[i]],val>0) 
  p1=ggplot(BaseError ,aes(val,color=sample))+
    geom_density()+scale_x_log10()+ggtitle(names(ListOfBaseErrorAllSamples)[i])+theme_bw();
  
  p2=ggplot(BaseError,aes(y=val,x=sample))+
    geom_violin(aes(fill=sample))+scale_y_log10()+ggtitle(names(ListOfBaseErrorAllSamples)[i])+theme_bw();
  
  
  ListOfSplitedErrors =split(ListOfBaseErrorAllSamples[[i]],as.character(ListOfBaseErrorAllSamples[[i]]$sample) )
  ListOfZeroRatios = lapply(ListOfSplitedErrors, FindZeroRatio)
  sapply(1:length(ListOfZeroRatios),function(i) ListOfZeroRatios[[i]]$sample <<- names(ListOfZeroRatios)[i]  )
  ZeroRatios = do.call(rbind,ListOfZeroRatios)
  
  
  p3 =ggplot(data=ZeroRatios, aes(y=Freq, x=sample, fill=variable)) +
    geom_bar(stat="identity",color='black',width=0.5)+theme_bw()
  
  gridExtra::grid.arrange(p1,p2,p3,nrow=1,ncol=3)
  
},simplify = F )
dev.off()



