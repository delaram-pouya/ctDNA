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




FitNormal <- lapply(ListOfSplitednonZeroErrorsLog, function(x) fitdistr(x, "normal") )
FitNormal <- lapply(FitNormal, function(x) x$estimate )
FitNormalMean <- lapply(FitNormal, function(x) x[['mean']])
FitNormalSD <- lapply(FitNormal, function(x) x[['sd']])


SAMPLE_SIZE = 1000
GeneratedErrorLog <- lapply(1:length(FitNormal),function(i) rnorm( SAMPLE_SIZE , FitNormalMean[[i]], FitNormalSD[[i]]))
GeneratedError <- lapply(1:length(FitNormal),function(i) 10^rnorm( SAMPLE_SIZE , FitNormalMean[[i]], FitNormalSD[[i]]))


Files <- list(FitNormal, FitNormalMean, FitNormalSD, GeneratedErrorLog, GeneratedError)
sapply(1:length(Files), 
       function(i) names(Files[[i]]) <<- names(ListOfSplitedErrors))

lapply(Files, function(x)lapply(x, head))




### Simulation Visualization

sapply(1:length(FitNormal), function(i) hist(x = rnorm( SAMPLE_SIZE , FitNormalMean[[i]], FitNormalSD[[i]]),
                                             main =names(FitNormal)[i], 
                                             xlab=names(FitNormal)[i]) )

sapply(1:length(FitNormal), function(i) hist(x = 10^(rnorm(SAMPLE_SIZE ,FitNormalMean[[i]], FitNormalSD[[i]])),
                                             main =names(FitNormal)[i], 
                                             xlab=names(FitNormal)[i] ) )  


ErrorSimulations <- sapply(1:length(FitNormal), function(i)  {
       data.frame( Sample=names(FitNormal)[i], 
                   Error=10^(rnorm(SAMPLE_SIZE ,FitNormalMean[[i]], FitNormalSD[[i]])) )
  }, simplify = F)
ErrorSimulations <- do.call(rbind, ErrorSimulations)

ggplot(data=ErrorSimulations, aes(x=Error)) +geom_density(alpha=0.3,aes(fill=Sample)) +theme_bw()+ scale_fill_discrete()


pdf('Error/Plots/EvalSimulate.pdf')
dev.off()



