
### importing data

ListOfSplitedErrorsMerged <- readRDS('error/Data/ListOfSplitedErrorsMerged.rds')
lapply(ListOfSplitedErrorsMerged, head)

ListOfSplitedErrorsMerged.df <- lapply(ListOfSplitedErrorsMerged, data.frame)
ListOfSplitedErrorsMerged.df <- lapply(ListOfSplitedErrorsMerged.df, function(x) {colnames(x)<-'error';x})
names(ListOfSplitedErrorsMerged.df) <- names(ListOfSplitedErrorsMerged)
samples <- names(ListOfSplitedErrorsMerged.df)



### prooving the need for log-transformation -> raw data is highly skewed

pdf('error/Plots/ditribution_raw_merged.pdf')
sapply(1:length(ListOfSplitedErrorsMerged.df), 
       function(i){
         p=ggplot(ListOfSplitedErrorsMerged.df[[i]], aes(x=error))+
           geom_density(fill='cyan',color='black',alpha=0.5)+theme_bw()+
           ggtitle(paste0(names(ListOfSplitedErrorsMerged.df)[i],' error distribution'))
         print(p)
         descdist(ListOfSplitedErrorsMerged.df[[i]]$error)
       })
dev.off()


## can't fit any distribution to the raw data(highly skewed) >> error
findIndivDistribution <- function(data){
  plot.legend <- c("normal",  "gamma", 'logistic') #"lognormal",
  fn <- fitdist(data, "norm")
  flg <- fitdist(data, "lnorm")
  fg <- fitdist(data, "gamma")
  fln <- fitdist(data, "logis")
  
  p1=denscomp(list(fn, fg, fln), legendtext = plot.legend)
  p2=qqcomp(list(fn, fg, fln), legendtext = plot.legend)
  p3=cdfcomp(list(fn, fg, fln), legendtext = plot.legend)
  p4=ppcomp(list(fn, fg, fln), legendtext = plot.legend)
  return(list(p1, p2, p3, p4))  
}



######## prooving that individual distribution can't work on the log-transformed data
##### checking if single gaussian fits better

pdf('Plots/errorDescdist_merged.pdf')
sapply(1:length(ListOfSplitednonZeroErrorsMergedLog), function(i){
  hist(ListOfSplitednonZeroErrorsMergedLog[[i]], main = names(ListOfSplitednonZeroErrorsMergedLog)[i] ,xlab = 'log10(non-zero error)') ;
  descdist(ListOfSplitednonZeroErrorsMergedLog[[i]], boot = 100)
},
simplify = F)

dev.off()


### checking CDF and QQ plots

samples <- c('blood', 'cfDNA', 'tissue')

sapply(1:length(samples), function(i){
  plot.legend <- c("normal", "lognormal", "gamma", 'logistic')
  fn <- fitdist(-ListOfSplitednonZeroErrorsMergedLog[[ samples[i] ]], "norm")
  flg <- fitdist(-ListOfSplitednonZeroErrorsMergedLog[[ samples[i] ]], "lnorm")
  fg <- fitdist(-ListOfSplitednonZeroErrorsMergedLog[[ samples[i] ]], "gamma")
  fln <- fitdist(-ListOfSplitednonZeroErrorsMergedLog[[ samples[i] ]], "logis")
  
  
  pdf(paste0('error/Plots/FitDist_merged/',samples[i],'_fitDist_plots.pdf'), width=15, height=15)
  par(mfrow = c(1, 1))
  denscomp(list(fn, flg, fg, fln), legendtext = plot.legend)
  qqcomp(list(fn, flg, fg, fln), legendtext = plot.legend)
  cdfcomp(list(fn, flg, fg, fln), legendtext = plot.legend)
  ppcomp(list(fn, flg, fg, fln), legendtext = plot.legend)
  dev.off()
  
}, simplify = F)




### Kolmogorov-Smirnov Test
#### all of them get rejected to fit to individual guassian distribution

SIZE= 10000
FitNormaltoMerged <- lapply(ListOfSplitednonZeroErrorsMergedLog, function(x) fitdistr(x, "normal") )

## cfDNA
ks.test(ListOfSplitednonZeroErrorsMergedLog[['cfDNA']],
        rnorm(SIZE, FitNormaltoMerged$cfDNA$estimate['mean'],FitNormaltoMerged$cfDNA$estimate['sd']) )

## Blood
ks.test(ListOfSplitednonZeroErrorsMergedLog[['blood']],
        rnorm(SIZE, FitNormaltoMerged$blood$estimate['mean'],FitNormaltoMerged$blood$estimate['sd']) )

## tissue
ks.test(ListOfSplitednonZeroErrorsMergedLog[['tissue']],
        rnorm(SIZE, FitNormaltoMerged$tissue$estimate['mean'],FitNormaltoMerged$tissue$estimate['sd']) )


