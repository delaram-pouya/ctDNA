### subtracing the small peak from blood and tissue from  
### the respective peak in the cfDNA sample  

set.seed(123)

source('error/Codes/4_findError_functions.R')
Initialize()
ListOfSplitednonZeroErrorsMergedLog <- readRDS('error/Data/ListOfSplitednonZeroErrorsMergedLog.rsd')
ListOfSplitedErrorsMerged <- readRDS('error/Data/ListOfSplitedErrorsMerged.rds')


findSmallPeakIndex <- function(mixEM_list){
  maxMu <- lapply(mixEM_list, function(x) max(x$mu[1], x$mu[2]))
  sapply(1:length(maxMu), function(i)ifelse(mixEM_list[[i]]$mu[1] ==maxMu[[i]], 1 , 2) )
}




#### load the fited mixture distributions:

normalmixEM_blood <- readRDS('error/Data/normalmixEM_blood.rds')
normalmixEM_tissue <- readRDS('error/Data/normalmixEM_tissue.rds')
normalmixEM_cfDNA <- readRDS('error/Data/normalmixEM_cfDNA.rds')


# finding the small peak(added PCR error) parameters

Names <- c('blood', 'tissue', 'cfDNA')
mixEM_list <- list(normalmixEM_blood, normalmixEM_tissue ,normalmixEM_cfDNA)

index = findSmallPeakIndex(mixEM_list)

mu_list <- sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$mu[ index[i]  ]  )
sigma_list <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$sigma[ index[i]  ]  )
pmix_list <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$lambda[ index[i]  ]  )

names(mu_list)= Names; names(sigma_list) = Names; names(pmix_list) = Names




## calculating the number of samples taken 

BLOOD_SIMULARED_READ_COUNT = 606160
TISSUE_SIMULARED_READ_COUNT = 35428
TOTAL_SIMULARED_READ_COUNT = BLOOD_SIMULARED_READ_COUNT + TISSUE_SIMULARED_READ_COUNT

BLOOD_READ_RATIO = BLOOD_SIMULARED_READ_COUNT/TOTAL_SIMULARED_READ_COUNT
TISSUE_READ_RATIO = TISSUE_SIMULARED_READ_COUNT/TOTAL_SIMULARED_READ_COUNT

## finding the non-zero ratio
nonZerolen <- unlist(lapply(ListOfSplitednonZeroErrorsMergedLog, length))
totalLen <- unlist(lapply(ListOfSplitedErrorsMerged, length))
NonZero_ratio <- nonZerolen/totalLen


# 1250103
SAMPLE_SIZE <- BLOOD_READ_RATIO * pmix_list[['blood']] * length(ListOfSplitedErrorsMerged[['blood']]) * NonZero_ratio[['blood']] +
  TISSUE_READ_RATIO * pmix_list[['tissue']] * length(ListOfSplitedErrorsMerged[['tissue']]) * NonZero_ratio[['tissue']]


### generating small peak samples for tissue and cfDNA and blood
pcr_error <- sapply(1:length(mixEM_list), function(i) 10**rnorm(SAMPLE_SIZE, mu_list[i], sigma_list[i]),simplify=F)
names(pcr_error) <- Names

added_tissue_error <- pcr_error[['cfDNA']] - pcr_error[['tissue']]
added_blood_error <- pcr_error[['cfDNA']] - pcr_error[['blood']]

summary(added_tissue_error)
summary(added_blood_error)

## few adjustments
added_blood_error[added_blood_error<0] = 0 
added_tissue_error[added_tissue_error<0] = 0
added_blood_error[added_blood_error>1] = 1 
added_tissue_error[added_tissue_error>1] = 1

summary(added_tissue_error)
summary(added_blood_error)


write.csv(x =added_blood_error, file = 'error/Data/bloodSmallPeak.csv',row.names = F, col.names = F, quote = F)
write.csv(x =added_tissue_error, file = 'error/Data/tissueSmallPeak.csv',row.names = F, col.names = F, quote = F)



lapply(mixEM_list, function(x) x$mu)
lapply(mixEM_list, function(x) x$sigma)
lapply(mixEM_list, function(x) x$lambda)
