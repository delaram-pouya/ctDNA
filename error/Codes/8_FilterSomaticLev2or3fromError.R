source('Error/Codes/findError_functions.R')
Initialize()



#### loading all pileup files
samples <- c('p1_blood','p1_cfDNA1','p1_cfDNA2','p1_tissue',
             'p2_blood','p2_cfDNA1','p2_cfDNA2','p2_tissue')


## HOMOZYGOUS 
#### finding homozygous positions for p1 & p2 based on blood sample
homo <- readRDS('Error/Data/homozygous_pos.rds')




######### TISSUE ##############
#### Somatic mutations

v1 <- fread('Error/Data/P1/VCFs/p1_tissue_snv_lev3.txt', header = T)
v2 <- fread('Error/Data/P2/VCFs/p2_tissue_snv_lev3.txt', header = T)
v <- list(v1,v2); names(v)=c('p1_tissue','p2_tissue')
somatics <- mclapply(v, function(x)paste0(x$chrom,sep='_',x$pos),mc.cores = detectCores()-2)


#### importing tissue pileup files 
tissueHomo <- readRDS('Error/Data/tisPileupHomoFilt.rds')


# remove somatic mutations
index_somatic <- sapply(1:length(tissueHomo), 
                        function(i) paste0(tissueHomo[[i]]$seqnames,sep='_',tissueHomo[[i]]$pos) %in% somatics[[i]],
                        simplify=F)

tissueHomoNoMut <- sapply( 1:length(tissueHomo), function(i)tissueHomo[[i]][ (!index_somatic[[i]]) ,],simplify = F )
tissueHomoNoMut <- AddAttributes(tissueHomoNoMut)
names(tissueHomoNoMut) <- c('p1_tissue','p2_tissue') 
tissueHomoNoMutError <- c(tissueHomoNoMut[[1]]$ALT.freq, tissueHomoNoMut[[2]]$ALT.freq)
tissueHomoNoMutNonZeroError <- tissueHomoNoMutError[tissueHomoNoMutError>0]
tissueHomoNoMutErrorDF <- data.frame(error=log10(tissueHomoNoMutNonZeroError), sample='Non-Somatics')
#saveRDS(tissueHomoNoMut, 'Data/tisPileupFinal_NoMutFilter.rds')



tissueMuts <- sapply( 1:length(tissueHomo), function(i)tissueHomo[[i]][ (index_somatic[[i]]) ,],simplify = F )
tissueMuts <- AddAttributes(tissueMuts)
names(tissueMuts) <- c('p1_tissue','p2_tissue') 
tissueMutsError <- c(tissueMuts[[1]]$ALT.freq, tissueMuts[[2]]$ALT.freq)
tissueMutsNonZeroError <- tissueMutsError[tissueMutsError>0]
tissueMutsErrorDF <- data.frame(error=log10(tissueMutsNonZeroError), sample='Somatics')



tissue_non_zero_homoFilt_all <- rbind(tissueHomoNoMutErrorDF, tissueMutsErrorDF)
p1=ggplot(tissue_non_zero_homoFilt_all, aes(error))+theme_bw()+ggtitle('Tissue')+geom_histogram(color='black',aes(fill=sample))
p2=ggplot(tissueMutsErrorDF, aes(error))+geom_histogram(color='black')+xlim(c(-5,0))+theme_bw()+ggtitle('Somatics')
p3=ggplot(tissueHomoNoMutErrorDF, aes(error))+geom_histogram(color='black')+xlim(c(-5,0))+theme_bw()+ggtitle('Non-Somatics')
grid.arrange(p1,p2,p3, nrow=1,ncol=3)







####  cfDNA error distributions ###

#### finding Somatic mutations

cfdna1VCF <- lapply(list.files('Error/Data/P1/VCFs/','lev2.txt$', full.names = T)[1:2], read.table, header=T)
cfdna2VCF <- lapply(list.files('Error/Data/P2/VCFs/','lev2.txt$', full.names = T)[1:2], read.table, header=T)
cfDNAvcf = c(cfdna1VCF,cfdna2VCF)
names(cfDNAvcf)=c('p1_cfDNA1','p1_cfDNA2','p2_cfDNA1','p2_cfDNA2')
cfDNAvcf <- lapply(cfDNAvcf, function(x){colnames(x)=c('chrom','pos'); return(x[2:nrow(x),]) })
somatics.cfDNA <- mclapply(cfDNAvcf, function(x)paste0(x$chrom,sep='_',x$pos),mc.cores = detectCores()-2)


#### loading cfDNA pileups 
# remove hetero positions
cfdnaHomo <- readRDS('Error/Data/cfdnaPileupHomoFilt.rds')


# remove somatic mutations
index_somatic <- sapply(1:length(cfdnaHomo), 
                        function(i) paste0(cfdnaHomo[[i]]$seqnames,sep='_',cfdnaHomo[[i]]$pos) %in% somatics.cfDNA[[i]],
                        simplify=F)

### Non-Somatic positions

cfdnaHomoNoMut <- sapply(1:length(cfdnaHomo), function(i) cfdnaHomo[[i]][ (!index_somatic[[i]]) ,],simplify = F )
cfdnaHomoNoMut <- AddAttributes(cfdnaHomoNoMut)
names(cfdnaHomoNoMut) <- c('p1_cfDNA1','p1_cfDNA2','p2_cfDNA1','p2_cfDNA2')
#saveRDS(cfdnaHomoNoMut, 'Data/cfDNAsPileupFinal_NoMutFilt.rds')
cfdnaHomoNoMutError <- c(cfdnaHomoNoMut[[1]]$ALT.freq, cfdnaHomoNoMut[[2]]$ALT.freq,
                         cfdnaHomoNoMut[[3]]$ALT.freq, cfdnaHomoNoMut[[4]]$ALT.freq)

cfdnaHomoNoMutError <- cfdnaHomoNoMutError[cfdnaHomoNoMutError>0]
cfdnaHomoNoMutErrorDF <- data.frame(error=log10(cfdnaHomoNoMutError), sample='non-Somatics')




### Somatic  positions

cfdnaHomoMuts <- sapply(1:length(cfdnaHomo), function(i) cfdnaHomo[[i]][ (index_somatic[[i]]) ,],simplify = F )
cfdnaHomoMuts <- AddAttributes(cfdnaHomoMuts)
names(cfdnaHomoMuts) <- c('p1_cfDNA1','p1_cfDNA2','p2_cfDNA1','p2_cfDNA2')
cfdnaHomoMutsError <- c(cfdnaHomoMuts[[1]]$ALT.freq, cfdnaHomoMuts[[2]]$ALT.freq,
                        cfdnaHomoMuts[[3]]$ALT.freq, cfdnaHomoMuts[[4]]$ALT.freq)

cfdnaHomoMutsError <- cfdnaHomoMutsError[cfdnaHomoMutsError>0]
cfdnaHomoMutsErrorDF <- data.frame(error=log10(cfdnaHomoMutsError), sample='Somatics')

cfDNA_non_zero_homoFilt_all <- rbind(cfdnaHomoNoMutErrorDF, cfdnaHomoMutsErrorDF)
p1=ggplot(cfDNA_non_zero_homoFilt_all, aes(error))+theme_bw()+ggtitle('cfDNA')+geom_histogram(color='black',aes(fill=sample))
p2=ggplot(cfdnaHomoMutsErrorDF, aes(error))+geom_histogram(color='black')+xlim(c(-5,0))+theme_bw()+ggtitle('Somatics')
p3=ggplot(cfdnaHomoNoMutErrorDF, aes(error))+geom_histogram(color='black')+xlim(c(-5,0))+theme_bw()+ggtitle('Non-Somatics')
grid.arrange(p1,p2,p3, nrow=1,ncol=3)




## checking number of removed somatic positions in each level 
### less than expected, some somcatic mutations in vcf file are not included in the pileup
#### removed by the filter

cfdnaHomolev2Filt <- readRDS('Error/Data/cfDNAsPileupFinal_lev2MutFilt.rds')
cfdnaHomolev3Filt <- readRDS('Error/Data/cfDNAsPileupFinal_lev3MutFilt.rds')
cfdnaHomoNumut <- readRDS('Error/Data/cfDNAsPileupFinal_NoMutFilt.rds')

unlist(lapply(cfdnaHomolev2Filt, nrow)) - unlist(lapply(cfdnaHomoNumut, nrow))
unlist(lapply(cfdnaHomolev3Filt, nrow)) - unlist(lapply(cfdnaHomoNumut, nrow))
lapply(somatics.cfDNA, length)




