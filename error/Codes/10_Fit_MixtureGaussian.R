source('error/Codes/4_findError_functions.R')
Initialize()
ListOfSplitednonZeroErrorsMergedLog <- readRDS('error/Data/ListOfSplitednonZeroErrorsMergedLog.rsd')


mergedNonzeroLog_blood <- ListOfSplitednonZeroErrorsMergedLog[['blood']]
#normalmixEM_blood = normalmixEM(mergedNonzeroLog_blood)
normalmixEM_blood <- readRDS('error/Data/normalmixEM_blood.rds')

pdf('error/Plots/normalmixEM_blood_plot.pdf')
plot(normalmixEM_blood,which=2)
lines(density(mergedNonzeroLog_blood), lty=2, lwd=2)
dev.off()


mergedNonzeroLog_tissue <- ListOfSplitednonZeroErrorsMergedLog[['tissue']]
#normalmixEM_tissue = normalmixEM(mergedNonzeroLog_tissue)
normalmixEM_tissue <- readRDS('error/Data/normalmixEM_tissue.rds')

pdf('error/Plots/normalmixEM_tissue_plot.pdf')
plot(normalmixEM_tissue,which=2)
lines(density(mergedNonzeroLog_tissue), lty=2, lwd=2)
dev.off()


mergedNonzeroLog_cfDNA <- ListOfSplitednonZeroErrorsMergedLog[['cfDNA']]
#normalmixEM_cfDNA = normalmixEM(mergedNonzeroLog_cfDNA)
normalmixEM_cfDNA <- readRDS('error/Data/normalmixEM_cfDNA.rds')

plot(normalmixEM_cfDNA ,which=2)
lines(density(mergedNonzeroLog_cfDNA), lty=2, lwd=2)





#### Goodness-of-fit test >> is this valid????

# CDF of mixture of two normals
pmnorm <- function(x, mu, sigma, pmix) {
  pmix[1]*pnorm(x,mu[1],sigma[1]) + (1-pmix[1])*pnorm(x,mu[2],sigma[2])
}

pnorm(normalmixEM_cfDNA$mu, 
      normalmixEM_cfDNA$sigma, 
      normalmixEM_cfDNA$lambda)


ks.cfDNA <- ks.test(mergedNonzeroLog_cfDNA, 
                pmnorm, 
                mu = normalmixEM_cfDNA$mu, 
                sigma = normalmixEM_cfDNA$sigma, 
                pmix = normalmixEM_cfDNA$lambda)

ks.blood <- ks.test(mergedNonzeroLog_blood, 
                    pmnorm, 
                    mu = normalmixEM_blood$mu, 
                    sigma = normalmixEM_blood$sigma, 
                    pmix = normalmixEM_blood$lambda)






#### ToDo:
# estimate cfdna + blood parameters by EM
# estimete tissue parameters (single guassian)
# Non-zero-SIZE = 10000
# 1000/total-size = non-zero ratio
# total-size * zero-ratio = number of zeros to add
# 10^value of samples taken from log + zeros >>> generated sample 
# take samples from the generated cfDNA sample - blood sample
# take samples from the generated cfDNA sample - tissue sample

## give these numbers to Ali 
## how many random positions should he choose???




ListOfSplitedErrorsMerged <- readRDS('Data/ListOfSplitedErrorsMerged.rds')
lapply(ListOfSplitedErrorsMerged, mean)

nonZerolen <- unlist(lapply(ListOfSplitednonZeroErrorsMergedLog, length))
totalLen <- unlist(lapply(ListOfSplitedErrorsMerged, length))
NonZero_ratio <- nonZerolen/totalLen
Zero_Ratio <- (1 - NonZero_ratio)

NUMBER_OF_SAMPLES_TO_GENERATE = 1000
NUMBER_OF_SAMPLES_TO_GENERATE * Zero_Ratio
NUMBER_OF_SAMPLES_TO_GENERATE * NonZero_ratio


