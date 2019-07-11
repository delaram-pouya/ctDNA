source('Error/Codes/findError_functions.R')
Initialize()

##########

tissueHomoNoMut <- readRDS('Error/Data/tisPileupFinal_lev2MutFilt.rds')
bloodPileupHomo <- readRDS('Error/Data/bloodfinalpileup.rds')
cfdnaHomoNoMut <- readRDS('Error/Data/cfDNAsPileupFinal_lev2MutFilt.rds')


TotalErrorLists <- rbind(
  data.frame(sample='tissue1',error=tissueHomoNoMut[[1]]$ALT.freq),
  data.frame(sample='tissue2',error=tissueHomoNoMut[[2]]$ALT.freq),
  
  data.frame(sample='blood1',error=bloodPileupHomo[[1]]$ALT.freq),
  data.frame(sample='blood2',error=bloodPileupHomo[[2]]$ALT.freq),
  
  data.frame(sample='p1_cfDNA1',error=cfdnaHomoNoMut[[1]]$ALT.freq),
  data.frame(sample='p1_cfDNA2',error=cfdnaHomoNoMut[[2]]$ALT.freq),
  data.frame(sample='p2_cfDNA1',error=cfdnaHomoNoMut[[3]]$ALT.freq),
  data.frame(sample='p2_cfDNA2',error=cfdnaHomoNoMut[[4]]$ALT.freq) )



################# Visualize

pdf('Analysis/errorAnalysis/errorAnalysisPlots/totalErrorAllSamples.pdf')
ggplot(TotalErrorLists,aes(error,color=sample))+geom_density()+
  ggtitle('error(F2+F3+f4)-homo-nonSNV')+theme_bw()
ggplot(TotalErrorLists,aes(error+1e-7,color=sample))+geom_density()+
  ggtitle('error(F2+F3+f4)-homo-nonSNV')+scale_x_log10()+xlab('error+1e-7(log)')+theme_bw()

ggplot(subset(TotalErrorLists,error>0),aes(error,color=sample))+geom_density()+
  ggtitle('error(F2+F3+f4)-homo-nonSNV-nonZero')+theme_bw()
ggplot(subset(TotalErrorLists,error>0),aes(error,color=sample))+geom_density()+
  ggtitle('error(F2+F3+f4)-homo-nonSNV-nonZero')+scale_x_log10()+xlab('error.log')+theme_bw()
dev.off()



# ZeroRatio all samples

colnames(TotalErrorLists)[2] <- 'val'
ListOfSplitedErrors =split(TotalErrorLists,as.character(TotalErrorLists$sample) )
ListOfZeroRatios = lapply(ListOfSplitedErrors, FindZeroRatio)
sapply(1:length(ListOfZeroRatios),function(i) ListOfZeroRatios[[i]]$sample <<- names(ListOfZeroRatios)[i]  )
ZeroRatios = do.call(rbind,ListOfZeroRatios)
ggplot(data=ZeroRatios, aes(y=Freq, x=sample, fill=variable)) +
  geom_bar(stat="identity",color='black',width=0.5)+theme_bw()



############## Error distribution parameter estimation
colnames(TotalErrorLists)[2] <- 'val'
TotalErrorListsNonZero <- subset(TotalErrorLists,val>0)

ListOfSplitedErrors  <- split(TotalErrorLists,as.character(TotalErrorLists$sample) )
ListOfSplitednonZeroErrors  <- split(TotalErrorListsNonZero,as.character(TotalErrorListsNonZero$sample) )
names(ListOfSplitednonZeroErrors) <- names(ListOfSplitedErrors)

ListOfSplitednonZeroErrorsLog <- lapply(ListOfSplitednonZeroErrors, function(x) log10(x$val))
names(ListOfSplitednonZeroErrorsLog) <- names(ListOfSplitedErrors)

pdf('Error/Plots/errorDescdist_lev2MutFilt.pdf')
sapply(1:length(ListOfSplitednonZeroErrorsLog), function(i){
  hist(ListOfSplitednonZeroErrorsLog[[i]], main = names(ListOfSplitednonZeroErrorsLog)[i] ,xlab = 'log10(error)') ;
  #descdist(ListOfSplitednonZeroErrorsLog[[i]])
  },
  simplify = F)

dev.off()



ErrorMeans <- lapply(ListOfSplitedErrors, function(x) mean(x$val))
NumberOfPositions <- lapply(ListOfSplitedErrors, nrow)
lapply(ListOfSplitednonZeroErrorsLog, nrow)
sapply(1:length(ErrorMeans), function(i) ErrorMeans[[i]]*NumberOfPositions[[i]])






##########################################################
## visualize error, based on each chromosome

tissueHomoNoMut <- readRDS('error/Data/tisPileupFinal_NoMutFilter.rds')
bloodPileupHomo <- readRDS('error/Data/bloodfinalpileup.rds')
cfdnaHomoNoMut <- readRDS('error/Data/cfDNAsPileupFinal_NoMutFilt.rds')


## adding a total plot(summation of all distributions)
cfdnaHomoNoMut <- lapply(cfdnaHomoNoMut, function(x){
  tmp = subset(x, seqnames != 'chrY') 
  tmp$seqnames <- 'total'
  return(rbind(x, tmp))
})




### make a wide range of colors to plot multiple classes
library(RColorBrewer)
n <- 50
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))



## same for error_blood_chr.pdf & error_cfDNA_chr.pdf

pdf('error/Plots/error_tissue_chr_2.pdf')
sapply(1:length(cfdnaHomoNoMut), 
       function(i) ggplot(cfdnaHomoNoMut[[i]], aes(x=ALT.freq+1e-7,color=seqnames))+geom_density()+theme_bw()+
         scale_x_log10()+ggtitle(paste0(names(cfdnaHomoNoMut),' heterozygous and somatics removed' ))+
         xlab('log scaled error(=alt-freq)')+scale_color_manual(values=col_vector)
       ,simplify = F)

dev.off()


