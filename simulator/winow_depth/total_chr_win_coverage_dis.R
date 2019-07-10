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
chr_coverage <- data.frame(chr_total)
names(chr_coverage) <- 'coverage'
chr_name <- 'total coverage'

print(paste0('analyzing ',chr_name))
pdf('simulator/winow_depth/Chromosome_counts/total_chr_distrib.pdf')

ggplot(chr_coverage,aes(x=coverage))+geom_histogram(color='black',fill='cyan',bins=45,alpha=0.6)+theme_bw()+ggtitle(chr_name)
ggplot(chr_coverage,aes(x=coverage))+geom_density(color='black',fill='pink',alpha=0.7)+theme_bw()+ggtitle(chr_name)
ggplot(chr_coverage,aes(x=coverage))+geom_density(color='black',fill='pink',alpha=0.5)+
  theme_bw()+ggtitle(paste0(chr_name,' log-transformed'))+scale_x_log10()


descdist(chr_coverage$coverage, discrete=T, boot=50)

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



poisson_parameters <- listOfDistr[[1]]$estimate
nbinomial_parameters <- listOfDistr[[2]]$estimate

## Goodness of fit test
## Chi-squared statistic 

## warning: approximation might be incorrect
pois <- stats::chisq.test(rpois(nrow(chr_coverage),  poisson_parameters), 
                          chr_coverage$coverage  )

nbinom <- stats::chisq.test(rnbinom(nrow(chr_coverage), size = nbinomial_parameters[1], mu=nbinomial_parameters[2] ), 
                            chr_coverage$coverage  )





