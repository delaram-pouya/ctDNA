## in this script we'll try to check the
## error distribution on the dataseries_1 data
## to check if the same pattern of distribution(bimodal) is consistent 

source('/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/ALN/depth_all/pileup_functions.R')
source('/media/pgdrive/sharif/cfdna/cfDNA/DataSeries2/ALN/depth_all/findError_functions.R')
Initialize()


AddAttributes_singleDF <- function(x){
  x <- preprop(x)
  x <- MakepileupFreq.eff(x)
  x <- Add.depth(x)
  x <- Add.max(x)
  return(x)}


setwd('/media/pgdrive/sharif/cfdna/cfDNA/Data_Series1/ALN/')

### importing variants
Somatic_lev2 <- list.files(path='/media/pgdrive/sharif/cfdna/cfDNA/Data_Series1/Analysis/Variants' , 
                           pattern = '*lev2.txt', full.names = T)
 
Somatic_lev2 <- lapply(Somatic_lev2, read.table, header=T)
Somatic_lev2 <- lapply(Somatic_lev2, function(x) paste0(x$chrom, sep='_', x$pos) )
names(Somatic_lev2) <- c('cfDNA1', 'cfDNA2', 'cfDNA3')

############
Somatic_lev3 <- list.files(path='/media/pgdrive/sharif/cfdna/cfDNA/Data_Series1/Analysis/Variants' , 
                           pattern = '*lev3.txt', full.names = T)

Somatic_lev3 <- lapply(Somatic_lev3, read.table, header=T)
Somatic_lev3 <- lapply(Somatic_lev3, function(x) paste0(x$chrom, sep='_', x$pos) )
names(Somatic_lev3) <- c('cfDNA1', 'cfDNA2', 'cfDNA3')


## blood 

P1_blood <- readRDS('Error/Data/D1_p1_bloodPileup.rds')
P1_blood <- data.table(P1_blood[[1]])
P1_blood2 <- AddAttributes_singleDF(P1_blood)
summary(P1_blood2$F1)
bloodHomoDF <- subset(P1_blood2,F1>0.8)
bloodHomoDF <- readRDS( 'Error/Data/D1_P1_blood_HomoFilt.rds')

### finding homozygous positions based on blood sample F1
bloodHomoPositions <- paste0(bloodHomoDF$seqnames,sep='_',bloodHomoDF$pos)
bloodHomoPositions <- readRDS('Error/Data/D1_Homo_positions.rds')

P1_Blood_Error <- bloodHomoDF$ALT.freq
P1_Blood_Error_nonZero <- P1_Blood_Error[P1_Blood_Error>0]
P1_Blood_Error_nonZero.log <- log10(P1_Blood_Error_nonZero)
ggplot(data.frame(P1_Blood_Error_nonZero.log), aes(P1_Blood_Error_nonZero.log))+geom_density()+theme_bw()





bloodHomoPositions <- readRDS('Error/Data/D1_Homo_positions.rds')

## cfDNA1 
cfDNA1_pileup <- readRDS('Error/Data/D1_p1_cfDNA1_Pileup.rds')
cfDNA1_pileup <- data.table(cfDNA1_pileup)
cfDNA1_pileup <- AddAttributes_singleDF(cfDNA1_pileup)
cfDNA1_homoIndex <- paste0(cfDNA1_pileup$seqnames ,sep='_',cfDNA1_pileup$pos) %in% bloodHomoPositions  
cfDNA1_HomoDF <- cfDNA1_pileup[cfDNA1_homoIndex,]    ## only positions present in blood are included
saveRDS(cfDNA1_HomoDF, 'Error/Data/D1_P1_cfDNA1_HomoFilt.rds')

## same as above
saveRDS(cfDNA2_HomoDF, 'Error/Data/D1_P1_cfDNA2_HomoFilt.rds') ## where are thsese files? on the master branch?
saveRDS(cfDNA3_HomoDF, 'Error/Data/D1_P1_cfDNA3_HomoFilt.rds')


### removing the somatic mutations
cfDNA1_somatic_lev2_Index <- !paste0(cfDNA1_HomoDF$seqnames ,sep='_',cfDNA1_HomoDF$pos) %in% Somatic_lev2[[1]]  
cfDNA1_HomoDF_nonSNV <- cfDNA1_HomoDF[cfDNA1_somatic_lev2_Index,]

cfDNA2_somatic_lev2_Index <- !paste0(cfDNA2_HomoDF$seqnames ,sep='_',cfDNA2_HomoDF$pos) %in% Somatic_lev2[[2]]  
cfDNA2_HomoDF_nonSNV <- cfDNA2_HomoDF[cfDNA2_somatic_lev2_Index,]

cfDNA3_somatic_lev2_Index <- !paste0(cfDNA3_HomoDF$seqnames ,sep='_',cfDNA3_HomoDF$pos) %in% Somatic_lev2[[3]]  
cfDNA3_HomoDF_nonSNV <- cfDNA3_HomoDF[cfDNA3_somatic_lev2_Index,]



SampleNames <- c('Blood', 'cfDNA1', 'cfDNA2', 'cfDNA3')
listOfSamples <- list(bloodHomoDF ,cfDNA1_HomoDF, cfDNA2_HomoDF, cfDNA3_HomoDF)

listOferrors <- lapply(listOfSamples, function(x) x$ALT.freq)
listOferrors <- lapply(listOferrors, function(x)log10(x[x>0]))
listOferrorsDF <- sapply(1:length(listOferrors), function(i) data.frame(error=listOferrors[[i]], sample=SampleNames[i]), simplify = F)
errorDF <- do.call(rbind, listOferrorsDF)


p1 = ggplot(errorDF, aes(x=error,color=sample))+geom_density()+theme_bw()+ggtitle('error-dis(homo)')

grid.arrange(p1, p2, nrow=1, ncol=2)



