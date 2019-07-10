#### In this script, we try to fit a distribution to the number of reads in 300-bp 
### window that is  sliding through the genome, 

Initialize <- function(){
  library("fitdistrplus")
  #library(actuar)
  library("grid")
  library("gridExtra")
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(qdapRegex)
}

Initialize()


#### coverage distribution 


## import p1_cfDNA1 data
PATH_1 = 'simulator/winow_depth/p1_cfDNA1_plots'
input_file_names_1  <- list.files(PATH_1,pattern = '*.txt',full.names = T,include.dirs = T)
input_file_names_1 <- input_file_names_1[-25]
input_files_1 <- lapply(input_file_names_1, read.table,header=F)
names(input_files_1) <- rm_between(input_file_names_1, "plots/p1_cfDNA1_", "_read", extract=TRUE)


## import p1_cfDNA2 data
PATH_2 = 'simulator/winow_depth/p1_cfDNA2_plots'
input_file_names_2  <- list.files(PATH_2,pattern = '*.txt',full.names = T,include.dirs = T)
input_file_names_2 <- input_file_names_2[-25]
input_files_2 <- lapply(input_file_names_2, read.table,header=F)
names(input_files_2) <- rm_between(input_file_names_2, "plots/p1_cfDNA2_", "_read", extract=TRUE)


## merging two files
input_files <- sapply(1:length(input_files_1), function(i) rbind(input_files_1[[i]], input_files_2[[i]]), simplify=F)
names(input_files) <- names(input_files_2)



input_files <- lapply(input_files, function(x) {colnames(x)<-c("window_num","read_freq");x})
chr_coverageList <- lapply(input_files, function(x){df=as.data.frame(x$read_freq);colnames(df)='coverage';df})
lapply(chr_coverageList, head)



## all chromosomes combined
read_freq_all <- lapply(input_files, function(x) x$read_freq)
names(read_freq_all) <- NULL
chr_total <- unlist(read_freq_all)
summary(chr_total)






distribution <- c('pois', 'nbinom')
distribution_name <- c("Poisson", "negative binomial")



FindTheDistribution <- function(i){
 
  chr_coverage = chr_coverageList[[i]]  
  chr_name = names(input_files)[i]
  
  print(paste0('analyzing ',chr_name))
  pdf(paste0('simulator/winow_depth/Chromosome_counts/',chr_name,'_distrib.pdf'))
   
  p1=ggplot(chr_coverage,aes(x=coverage))+geom_histogram(color='black',fill='cyan',bins=45,alpha=0.6)+theme_bw()+ggtitle(chr_name)
  p2=ggplot(chr_coverage,aes(x=coverage))+geom_density(color='black',fill='pink',alpha=0.7)+theme_bw()+ggtitle(chr_name)
  p3=ggplot(chr_coverage,aes(x=coverage))+geom_density(color='black',fill='pink',alpha=0.5)+
    theme_bw()+ggtitle(paste0(chr_name,' log-transformed'))+scale_x_log10()
  
  print(p1)
  print(p2)
  print(p3)
  
  descdist(chr_coverage$coverage, discrete=T, boot=500)
  
  listOfDistr <- sapply(1:length(distribution), function(i){
    ## fitting negative binomial/ poisson
    distib <- fitdist(data = chr_coverage$coverage, distr = distribution[i], discrete = T)
    print(plot(distib))
    grid.text(distribution_name[i] ,x = 0.36, y=0.8 ,gp=gpar(fontsize=9, col="black"),just = "center")
    grid.text(distribution_name[i] ,x = 0.85, y=0.8 ,gp=gpar(fontsize=9, col="black"),just = "center")
    print(summary(distib))
    return(distib)}
    ,simplify = F)
  
  ## compare CDF of two distributions    
  cdfcomp(listOfDistr)
  dev.off()
  
  ## Goodness of fit test
  ## Chi-squared statistic 
  
  SAMPLE_SIZE = 1000
  pois <- stats::chisq.test(rpois(SAMPLE_SIZE,  listOfDistr[[1]]$estimate[1]), 
                            chr_coverage$coverage[round(runif(SAMPLE_SIZE,1,length(chr_coverage$coverage)))]  )
  
  nbinom <- stats::chisq.test(rnbinom(SAMPLE_SIZE, size = listOfDistr[[2]]$estimate[1], mu=listOfDistr[[2]]$estimate[2] ), 
                    chr_coverage$coverage[round(runif(SAMPLE_SIZE,1,length(chr_coverage$coverage)))]  )
  
  return(list(listOfDistr, pois,nbinom))
}




## fitting distribution, estimating parameters and checking with GOF test
fitResults <- sapply(1:length(chr_coverageList), function(i)FindTheDistribution(i), simplify = F)
saveRDS(fitResults, 'simulator/winow_depth/chrDisFitResults.rds')
names(fitResults) <- names(input_files)



## comparing the p values for nbinom and poise distributions
pois_pvalues <- sapply(1:length(input_files), function(i)fitResults[[i]][[2]]$p.value)
nbinom_pvalues <- sapply(1:length(input_files), function(i)fitResults[[i]][[3]]$p.value)

chi_test <- data.frame(pois=pois_pvalues,nbinom=nbinom_pvalues)
rownames(chi_test) <- names(input_files)
chi_test$dis <- ifelse( (chi_test$nbinom<0.05) & (chi_test$pois)>0.05, 'pois', 'nbinom')
write.csv(chi_test, 'simulator/winow_depth/chi_test.csv')






####### Sample generation

generateSample <- function(i){

    if(chi_test$dis[i]=='nbinom'){
    sample = rnbinom(nrow(chr_coverageList[[i]]),
                     size=fitResults[[i]][[1]][[2]]$estimate[1], 
                     mu=fitResults[[i]][[1]][[2]]$estimate[2] )
    
  }else sample = rpois(nrow(chr_coverageList[[i]]),lambda = fitResults[[i]][[1]][[1]]$estimate[1] )

  return(sample)
}



getParameters <- function(i){
  if(chi_test$dis[i]=='nbinom'){
    parameter = fitResults[[i]][[1]][[2]]$estimate
  
    }else parameter = fitResults[[i]][[1]][[1]]$estimate
  
  return(parameter)
}



## generating samples for each chromosome
samples <- sapply(1:length(chr_coverageList), function(i) generateSample(i), simplify = F)

## writing samples to csv files
sapply(1:length(chr_coverageList), 
       function(i) 
         write.csv(samples[[i]],paste0('simulator/winow_depth/chr_counts_corrected_sampling/',names(input_files)[i],'_samples.csv')))


parameters <- sapply(1:length(chr_coverageList), function(i) getParameters(i), simplify = F)
names(parameters) <- names(input_files)
saveRDS(parameters, 'simulator/winow_depth/chr_counts_corrected_sampling/parameters.rds')
  

