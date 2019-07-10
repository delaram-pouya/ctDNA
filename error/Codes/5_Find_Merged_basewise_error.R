#### in this script, we'll try to compare the base-wise error distribution of the merged samples
## input: base-wise error tables 
## output: base-wise error distribution plots for merged sample types (all-cfDNAs, all-tissue, all-bloods)


source('Error/Codes/4_findError_functions.R')
Initialize()

# c('p1_cfDNA1','p1_cfDNA2','p2_cfDNA1','p2_cfDNA2')

## Base-wise error distributions 

## cfDNA
cfDNAMat1 <- readRDS('Error/Data/ErrorMatrix/cfDNAMat1.rds')
cfDNAMat2 <- readRDS('Error/Data/ErrorMatrix/cfDNAMat2.rds')
cfDNAMat3 <- readRDS('Error/Data/ErrorMatrix/cfDNAMat3.rds')
cfDNAMat4 <- readRDS('Error/Data/ErrorMatrix/cfDNAMat4.rds')
listOfcfDNAs <- list(cfDNAMat1, cfDNAMat2, cfDNAMat3, cfDNAMat4)
merged_cfDNA_muts <- lapply(listOfcfDNAs, makeErrorTable )
lapply(merged_cfDNA_muts, head)
merged_cfDNA_muts <- do.call(rbind, merged_cfDNA_muts)
dim(merged_cfDNA_muts)
makePlot(merged_cfDNA_muts, 'cfDNA')

## blood
bloodMat1 <- readRDS('Error/Data/ErrorMatrix/bloodMat1.rds')
bloodMat2 <- readRDS('Error/Data/ErrorMatrix/bloodMat2.rds')
listOfbloods <- list(bloodMat1, bloodMat2)
merged_blood_muts <- lapply(listOfbloods, makeErrorTable )
lapply(merged_blood_muts, head)
merged_blood_muts <- do.call(rbind, merged_blood_muts)
dim(merged_blood_muts)
makePlot(merged_blood_muts, 'blood')

## tissue
tissueMat1 <- readRDS('Error/Data/ErrorMatrix/tissueMat1.rds')
tissueMat2 <- readRDS('Error/Data/ErrorMatrix/tissueMat2.rds')
listOftissue <- list(tissueMat1, tissueMat2)
merged_tissue_muts <- lapply(listOftissue, makeErrorTable )
lapply(merged_tissue_muts, head)
merged_tissue_muts <- do.call(rbind, merged_tissue_muts)
dim(merged_tissue_muts)
makePlot(merged_tissue_muts, 'tissue')




############ making a list of error dataframes
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
                                     data.frame(val=baseWiseErrorTissue1[[i]],sample='tissue'),
                                     data.frame(val=baseWiseErrorTissue2[[i]],sample='tissue'),
                                     
                                     data.frame(val=baseWiseErrorBlood1[[i]],sample='blood'),
                                     data.frame(val=baseWiseErrorBlood2[[i]],sample='blood'), 
                                     
                                     data.frame(val=baseWiseError_cfDNA1[[i]],sample='cfDNA'), 
                                     data.frame(val=baseWiseError_cfDNA2[[i]],sample='cfDNA'), 
                                     data.frame(val=baseWiseError_cfDNA3[[i]],sample='cfDNA'), 
                                     data.frame(val=baseWiseError_cfDNA4[[i]],sample='cfDNA') 
                                   )}
                                   ,simplify = F)

names(ListOfBaseErrorAllSamples) = names(baseWiseErrorBlood1)


######## compare different samples for each error  
#### including Zeros

pdf('error/Plots/ListOfBaseErrorMergedSamplesZeroIncluded.pdf',width=14,height=7)
sapply(1:length(ListOfBaseErrorAllSamples),function(i){
  
  p1=ggplot(ListOfBaseErrorAllSamples[[i]],aes(val+1e-7,color=sample))+
    geom_density()+scale_x_log10()+ggtitle(names(ListOfBaseErrorAllSamples)[i])+theme_bw();
  
  p2=ggplot(ListOfBaseErrorAllSamples[[i]],aes(y=val+1e-7,x=sample))+
    geom_violin(aes(fill=sample))+scale_y_log10()+ggtitle(names(ListOfBaseErrorAllSamples)[i])+theme_bw();
  
  gridExtra::grid.arrange(p1,p2,nrow=1,ncol=2)
},simplify = F )
dev.off()


######## Non-Zero errors

pdf('Error/Plots/ListOfBaseErrorMergedSamplesNonZero.pdf',width=23,height = 6)
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







