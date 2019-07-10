source('error/Codes/4_findError_functions.R')
Initialize()

##########

tissueHomoNoMut <- readRDS('error/Data/tisPileupFinal_NoMutFilter.rds')
bloodPileupHomo <- readRDS('error/Data/bloodfinalpileup.rds')
cfdnaHomoNoMut <- readRDS('error/Data/cfDNAsPileupFinal_NoMutFilt.rds')


TotalErrorLists <- rbind(
  data.frame(sample='tissue1',error=tissueHomoNoMut[[1]]$ALT.freq),
  data.frame(sample='tissue2',error=tissueHomoNoMut[[2]]$ALT.freq),
  
  data.frame(sample='blood1',error=bloodPileupHomo[[1]]$ALT.freq),
  data.frame(sample='blood2',error=bloodPileupHomo[[2]]$ALT.freq),
  
  data.frame(sample='p1_cfDNA1',error=cfdnaHomoNoMut[[1]]$ALT.freq),
  data.frame(sample='p1_cfDNA2',error=cfdnaHomoNoMut[[2]]$ALT.freq),
  data.frame(sample='p2_cfDNA1',error=cfdnaHomoNoMut[[3]]$ALT.freq),
  data.frame(sample='p2_cfDNA2',error=cfdnaHomoNoMut[[4]]$ALT.freq) )


############## Error distribution parameter estimation
colnames(TotalErrorLists)[2] <- 'val'
TotalErrorListsNonZero <- subset(TotalErrorLists,val>0)

ListOfSplitedErrors  <- split(TotalErrorLists,as.character(TotalErrorLists$sample) )
ListOfSplitednonZeroErrors  <- split(TotalErrorListsNonZero,as.character(TotalErrorListsNonZero$sample) )
names(ListOfSplitednonZeroErrors) <- names(ListOfSplitedErrors)

ListOfSplitednonZeroErrorsLog <- lapply(ListOfSplitednonZeroErrors, function(x) log10(x$val))
names(ListOfSplitednonZeroErrorsLog) <- names(ListOfSplitedErrors)


pdf('Error/Plots/errorDescdist_NoMutFilter.pdf')
sapply(1:length(ListOfSplitednonZeroErrorsLog), function(i){
  hist(ListOfSplitednonZeroErrorsLog[[i]], main = names(ListOfSplitednonZeroErrorsLog)[i] ,xlab = 'log10(error)') ;
  descdist(ListOfSplitednonZeroErrorsLog[[i]])
},
simplify = F)
dev.off()


##### Merging samples

ListOfSplitedErrorsMerged <- list(c(ListOfSplitedErrors[[1]]$val, ListOfSplitedErrors[[2]]$val),
                                  c(ListOfSplitedErrors[[3]]$val, ListOfSplitedErrors[[4]]$val, ListOfSplitedErrors[[5]]$val, ListOfSplitedErrors[[6]]$val),
                                  c(ListOfSplitedErrors[[7]]$val, ListOfSplitedErrors[[8]]$val) )


names(ListOfSplitedErrorsMerged) <- c('blood', 'cfDNA', 'tissue')                                  
#saveRDS(ListOfSplitedErrorsMerged, 'Error/Data/ListOfSplitedErrorsMerged.rds')


###### Merged section

ListOfSplitedErrorsMerged <- readRDS('error/Data/ListOfSplitedErrorsMerged.rds')
ListOfSplitedErrorsMerged.df <- rbind(data.frame(sample='blood',error=ListOfSplitedErrorsMerged[['blood']]),
                                      data.frame(sample='cfDNA',error=ListOfSplitedErrorsMerged[['cfDNA']]),
                                      data.frame(sample='tissue',error=ListOfSplitedErrorsMerged[['tissue']])
                                      )
### zero-included hisogram
ggplot(ListOfSplitedErrorsMerged.df, aes(x=error))+
  geom_histogram(aes(fill=sample),color='black',alpha=0.6, bins=45)+scale_colour_manual(values=c("red","green","blue"))+
  theme_bw()+ggtitle('Error histogram-zero included')


### non-zero hisogram
ListOfSplitedErrorsMerged_nonzero.df <- subset(ListOfSplitedErrorsMerged.df,error>0)
ggplot(ListOfSplitedErrorsMerged_nonzero.df, aes(x=error))+
  geom_histogram(aes(fill=sample),color='black', alpha=0.6, bins=45)+scale_colour_manual(values=c("red","green","blue"))+
  theme_bw()+ggtitle('Error histogram-non zero')

## non-zero density
ggplot(ListOfSplitedErrorsMerged_nonzero.df, aes(x=error))+
  geom_density(aes(fill=sample), alpha=0.6)+scale_colour_manual(values=c("red","green","blue"))+
  theme_bw()+ggtitle('Error density non zero')



lapply(ListOfSplitedErrorsMerged, function(x) sum(x==0)/length(x) )
ListOfSplitednonZeroErrorsMerged <- lapply(ListOfSplitedErrorsMerged, function(x) x[x>0])
names(ListOfSplitednonZeroErrorsMerged) <- names(ListOfSplitedErrorsMerged) 


ListOfSplitednonZeroErrorsMergedLog <- lapply(ListOfSplitednonZeroErrorsMerged, function(x) log10(x))
names(ListOfSplitednonZeroErrorsMergedLog) <- names(ListOfSplitednonZeroErrorsMerged)
saveRDS(ListOfSplitednonZeroErrorsMergedLog, 'error/Data/ListOfSplitednonZeroErrorsMergedLog.rsd')



ErrorMeans <- lapply(ListOfSplitedErrorsMerged, function(x) mean(x))

lapply(ListOfSplitedErrorsMerged, function(x) summary(x))
lapply(ListOfSplitednonZeroErrorsMerged, function(x) summary(x))

MergedErrorNames <- names(ListOfSplitednonZeroErrorsMergedLog)


mergedLognonZeroErrorDF <-rbind( data.frame(error=ListOfSplitednonZeroErrorsMergedLog[[1]],
                                      sample= MergedErrorNames[1]), 
                                 data.frame(error=ListOfSplitednonZeroErrorsMergedLog[[2]],
                                            sample= MergedErrorNames[2]),
                                 data.frame(error=ListOfSplitednonZeroErrorsMergedLog[[3]],
                                            sample= MergedErrorNames[3])
                                 )

saveRDS(mergedLognonZeroErrorDF, 'error/Data/mergedLognonZeroErrorDF.rds')
head(mergedLognonZeroErrorDF)
dim(mergedLognonZeroErrorDF)
p1 = ggplot(mergedLognonZeroErrorDF, aes(x=error,color=sample))+geom_density()+theme_bw()+ggtitle('Error density non-zero (log scaled)')


zeroRationMerged <- melt(data.frame(lapply(ListOfSplitedErrorsMerged, function(x) sum(x==0)/length(x) )))
p2 = ggplot(zeroRationMerged, aes(y=value,x=variable))+
  geom_bar(stat="identity",width=0.5 ,aes(fill=variable),color='black')+
  theme_bw()+xlab('sample')+ylab('zero-Ratio')+ggtitle('zero-error-rate frequency')



pdf('Error/Plots/mergedErrorDis.pdf', width=15)
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()



