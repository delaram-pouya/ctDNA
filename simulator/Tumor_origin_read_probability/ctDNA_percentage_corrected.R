Initialize <- function(){
  library("fitdistrplus")
  #library(actuar)
  library("grid")
  library("gridExtra")
  library(dplyr)
  library(ggplot2)
  library(stringr)}

Initialize()


## Import data
pileupVAF <- readRDS('simulator/Tumor_origin_read_probability/pileUpVAFall.rds')
vcf <- readRDS('simulator/Tumor_origin_read_probability/allVCFs.rds')

lapply(pileupVAF, head)
lapply(vcf, head)

## merging data
merged_VAFs <- c(pileupVAF[['p1_cfDNA1']]$VAF,
                 pileupVAF[['p1_cfDNA2']]$VAF, 
                 pileupVAF[['p2_cfDNA1']]$VAF,
                 pileupVAF[['p2_cfDNA2']]$VAF)


merged_VAFs <- data.frame(merged_VAFs[merged_VAFs>0])
colnames(merged_VAFs) <- 'VAF'



pdf('simulator/Tumor_origin_read_probability/ctDNA_perc_lev2_corrected.pdf')
ggplot(data=merged_VAFs, aes(VAF)) + 
  geom_histogram(col="blue", fill="#56B4E9", alpha = .2 ,bins=20) + 
  labs(title="tumor-origin reads distribution") +
  labs(x=" %confirming alterations", y="frequency") +theme_bw()


ggplot(data=merged_VAFs, aes(VAF)) + 
  geom_density(col="blue", fill="#56B4E9", alpha = .2) + 
  labs(title="tumor-origin reads distribution") +
  labs(x=" %confirming alterations", y="frequency") +theme_bw()



x = merged_VAFs$VAF
x_scaled <- (x - min(x) + 0.0001) / (max(x) - min(x) + 0.0002)
descdist(x_scaled, discrete=F, boot=500)

fitParameters <- .FitDistributions(x_scaled)
dev.off()

saveRDS(fitParameters, 'simulator/Tumor_origin_read_probability/fitParameters_lev2_corrected.rds')




.FitDistributions <- function(ctData_perc){
  
  ## fit beta distribution
  beta_dis <- fitdist(data = ctData_perc, distr = "beta", discrete = F, method = "mle")
  plot(beta_dis)
  grid.text("Beta distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")
  
  ## fit normal distribution
  normal_dist <- fitdist(ctData_perc, "norm")
  plot(normal_dist)
  grid.text("Normal distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")
  
  ## fit gamma distribution
  gamma_dis <- fitdist(data = ctData_perc, distr = "gamma", discrete = F, method = "mle")
  plot(gamma_dis)
  grid.text("gamma distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")
  
  ## fit lognormal distribution
  lnorm_dis <- fitdist(data = ctData_perc, distr = "lnorm", discrete = F, method = "mle")
  plot(lnorm_dis)
  grid.text("lognorm distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")
  
  ## fit weibull distribution
  weibull_dis <- fitdist(data = ctData_perc, distr = "weibull", discrete = F, method = "mle")
  plot(weibull_dis)
  grid.text("weibull distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")
  
  
  legends <- c("beta", "norm" , "gamma", "lnorm", "weibull") 
  listOfparam <- list(beta_dis, normal_dist, gamma_dis, lnorm_dis, weibull_dis)
  names(listOfparam) <- legends
  
  denscomp(listOfparam, legendtext = legends)
  cdfcomp(listOfparam, legendtext = legends)
  qqcomp(listOfparam, legendtext = legends)
  ppcomp(listOfparam, legendtext = legends)
  par(mfrow=c(1,1)) 
  
  return(listOfparam)
}




# goodness-of-it statistics:
## choosing gamma as the best fit

gof = gofstat(fitParameters, fitnames = c("beta",  "norm" , "gamma", "lnorm", "weibull"))
gof$discrete

SAMPLE_SIZE = 500
ks.test(rbeta(SAMPLE_SIZE, fitParameters$beta$estimate[1],fitParameters$beta$estimate[2]),
        merged_VAFs$VAF[round(runif(SAMPLE_SIZE,1,length(merged_VAFs$VAF)))] )

ks.test(rnorm(SAMPLE_SIZE, fitParameters$norm$estimate[1],fitParameters$norm$estimate[2]),
        merged_VAFs$VAF[round(runif(SAMPLE_SIZE,1,length(merged_VAFs$VAF)))] )

ks.test(rgamma(SAMPLE_SIZE, fitParameters$gamma$estimate[1],fitParameters$gamma$estimate[2]),
        merged_VAFs$VAF[round(runif(SAMPLE_SIZE,1,length(merged_VAFs$VAF)))] )

ks.test(rlnorm(SAMPLE_SIZE, fitParameters$lnorm$estimate[1],fitParameters$lnorm$estimate[2]),
        merged_VAFs$VAF[round(runif(SAMPLE_SIZE,1,length(merged_VAFs$VAF)))] )

ks.test(rweibull(SAMPLE_SIZE, fitParameters$weibull$estimate[1],fitParameters$weibull$estimate[2]),
        merged_VAFs$VAF[round(runif(SAMPLE_SIZE,1,length(merged_VAFs$VAF)))] )




####### Sampling
num.reads = 86612989
Win.size = 300
win.count = 1000 

summary(rgamma(1000, shape = fitParameters$gamma$estimate[1] ,rate = fitParameters$gamma$estimate[2]))






###  read count distribution 
nbinom_sample =  rnbinom(win.count , size = 0.215427, mu =  156.116500)
summary(nbinom_sample)
hist(nbinom_sample)
write.csv(nbinom_sample, file = "C:/Users/Delaram/Desktop/coverage_sample.csv")



