## evaluate and compare small and large peaks in different samples
## using the KS test

set.seed(123)

source('error/Codes/4_findError_functions.R')
Initialize()
ListOfSplitednonZeroErrorsMergedLog <- readRDS('error/Data/ListOfSplitednonZeroErrorsMergedLog.rsd')
ListOfSplitedErrorsMerged <- readRDS('error/Data/ListOfSplitedErrorsMerged.rds')


findLargePeakIndex <- function(mixEM_list){
  minMu <- lapply(mixEM_list, function(x) min(x$mu[1], x$mu[2]))
  sapply(1:length(minMu), function(i)ifelse(mixEM_list[[i]]$mu[1] ==minMu[[i]], 1 , 2) )
}


findSmallPeakIndex <- function(mixEM_list){
  maxMu <- lapply(mixEM_list, function(x) max(x$mu[1], x$mu[2]))
  sapply(1:length(maxMu), function(i)ifelse(mixEM_list[[i]]$mu[1] ==maxMu[[i]], 1 , 2) )
}


#### load the fited mixture distributions:

normalmixEM_blood <- readRDS('error/Data/normalmixEM_blood.rds')
normalmixEM_tissue <- readRDS('error/Data/normalmixEM_tissue.rds')
normalmixEM_cfDNA <- readRDS('error/Data/normalmixEM_cfDNA.rds')


# finding the large peak(added PCR error) parameters

Names <- c('blood', 'tissue', 'cfDNA')
mixEM_list <- list(normalmixEM_blood, normalmixEM_tissue ,normalmixEM_cfDNA)





############## large peak 

index = findLargePeakIndex(mixEM_list)

mu_list <- sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$mu[ index[i]  ]  )
sigma_list <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$sigma[ index[i]  ]  )
pmix_list <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$lambda[ index[i]  ]  )

names(mu_list)= Names; names(sigma_list) = Names; names(pmix_list) = Names

SAMPLE_SIZE = 100
ks.test(rnorm(SAMPLE_SIZE,mu_list['blood'],sigma_list['blood']),
        rnorm(SAMPLE_SIZE,mu_list['tissue'],sigma_list['tissue']),alternative = 'two.sided')

ks.test(rnorm(SAMPLE_SIZE,mu_list['blood'],sigma_list['blood']),
        rnorm(SAMPLE_SIZE,mu_list['cfDNA'],sigma_list['cfDNA']),alternative = 'two.sided')

ks.test(rnorm(SAMPLE_SIZE,mu_list['tissue'],sigma_list['tissue']),
        rnorm(SAMPLE_SIZE,mu_list['cfDNA'],sigma_list['cfDNA']),alternative = 'two.sided')




############## small peak

index = findSmallPeakIndex(mixEM_list)

mu_list <- sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$mu[ index[i]  ]  )
sigma_list <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$sigma[ index[i]  ]  )
pmix_list <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$lambda[ index[i]  ]  )

names(mu_list)= Names; names(sigma_list) = Names; names(pmix_list) = Names

SAMPLE_SIZE = 100
ks.test(rnorm(SAMPLE_SIZE,mu_list['blood'],sigma_list['blood']),
        rnorm(SAMPLE_SIZE,mu_list['tissue'],sigma_list['tissue']),alternative = 'two.sided')

ks.test(rnorm(SAMPLE_SIZE,mu_list['blood'],sigma_list['blood']),
        rnorm(SAMPLE_SIZE,mu_list['cfDNA'],sigma_list['cfDNA']),alternative = 'two.sided')

ks.test(rnorm(SAMPLE_SIZE,mu_list['tissue'],sigma_list['tissue']),
        rnorm(SAMPLE_SIZE,mu_list['cfDNA'],sigma_list['cfDNA']),alternative = 'two.sided')




#######################################



############## large peak 

large_index = findLargePeakIndex(mixEM_list)

mu_large_peak <- sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$mu[ large_index[i]  ]  )
sigma_large_peak <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$sigma[ large_index[i]  ]  )
pmix_large_peak <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$lambda[ large_index[i]  ]  )




############## small peak

small_index = findSmallPeakIndex(mixEM_list)

mu_small_peak <- sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$mu[ small_index[i]  ]  )
sigma_small_peak <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$sigma[ small_index[i]  ]  )
pmix_small_peak <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$lambda[ small_index[i]  ]  )


names(mu_large_peak)= Names; names(sigma_large_peak) = Names; names(pmix_large_peak) = Names
names(mu_small_peak)= Names; names(sigma_small_peak) = Names; names(pmix_small_peak) = Names




n = 100000





