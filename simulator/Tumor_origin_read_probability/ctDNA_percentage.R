Initialize <- function(){
  library("fitdistrplus")
  #library(actuar)
  library("grid")
  library("gridExtra")
  library(dplyr)
  library(ggplot2)
  library(stringr)}

Initialize()


#### import data
ctData_perc_lev3 = read.csv("simulator/Tumor_origin_read_probability/p1_cfDNA1_lev3_ctdna_percent.txt", header = F)
ctData_perc_lev2 = read.csv("simulator/Tumor_origin_read_probability/p1_cfDNA1_lev2_ctdna_percent.txt", header = F)

ctData_perc <- ctData_perc_lev2
colnames(ctData_perc) <- 'tumor'

pdf('Simulator/Tumor_origin_read_probability/ctDNA_perc_lev2.pdf')
ggplot(data=ctData_perc, aes(ctData_perc$tumor)) + 
  geom_histogram(col="blue", fill="#56B4E9", alpha = .2 ,bins=20) + 
  labs(title="tumor-origin reads distribution") +
  labs(x=" %confirming alterations", y="frequency") +theme_bw()

ggplot(data=ctData_perc, aes(ctData_perc$tumor)) + 
  geom_density(col="blue", fill="#56B4E9", alpha = .2 ,bins=8) + 
  labs(title="tumor-origin reads distribution") +
  labs(x=" %confirming alterations", y="frequency") +theme_bw()

summary(ctData_perc)
descdist(ctData_perc$tumor/100, discrete=F, boot=500)
fitParameters <- .FitDistributions(ctData_perc)
dev.off()

saveRDS(fitParameters, 'Simulator/Tumor_origin_read_probability/fitParameters_lev2.rds')


.FitDistributions <- function(ctData_perc){
  
  ## fit beta distribution
  beta_dis <- fitdist(data = ctData_perc$tumor/100, distr = "beta", discrete = F, method = "mle")
  plot(beta_dis)
  grid.text("Beta distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")
  
  ## fit normal distribution
  normal_dist <- fitdist(ctData_perc$tumor/100, "norm")
  plot(normal_dist)
  grid.text("Normal distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")

  ## fit gamma distribution
  gamma_dis <- fitdist(data = ctData_perc$tumor/100, distr = "gamma", discrete = F, method = "mle")
  plot(gamma_dis)
  grid.text("gamma distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")
  
  ## fit lognormal distribution
  lnorm_dis <- fitdist(data = ctData_perc$tumor/100, distr = "lnorm", discrete = F, method = "mle")
  plot(lnorm_dis)
  grid.text("lognorm distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")

  ## fit weibull distribution
  weibull_dis <- fitdist(data = ctData_perc$tumor/100, distr = "weibull", discrete = F, method = "mle")
  plot(weibull_dis)
  grid.text("weibull distribution",x = 0.36,y= 0.85 ,gp=gpar(fontsize=9, col="black"),just = "center")

  
  legends <- c("beta", "norm" , "gamma", "lnorm", "weibull")
  listOfparam <- list(beta_dis, normal_dist, gamma_dis, lnorm_dis, weibull_dis)
  names(listOfparam) <- legends
  
  denscomp(listOfparam, legendtext = legends)
  cdfcomp(listOfparam, legendtext = legends)
  qqcomp(listOfparam, legendtext = legends)
  ppcomp(listOfparam, legendtext = legends)
  
  return(listOfparam)
}


# goodness-of-it statistics:
gof = gofstat(fitParameters, fitnames = c("beta", "norm" , "gamma", "lnorm", "weibull"))
gof



####### Sampling
num.reads = 86612989
Win.size = 300
win.count = 1000 ## ?


# Parameters : 
#   estimate Std. Error
# shape1 0.6745956 0.05415342
# shape2 1.4104338 0.13201152

beta_sample = rbeta(win.count, fitParameters$beta$estimate[1], fitParameters$beta$estimate[2])
summary(beta_sample)
hist(beta_sample)
#write.csv(beta_sample, file = "C:/Users/Delaram/Desktop/ctPercentage_sample.csv")


# Parameters : 
#         estimate   Std. Error 
# size   0.215427 0.0002586454
# mu   156.116500 0.3670577503
nbinom_sample =  rnbinom(win.count , size = 0.215427, mu =  156.116500)
summary(nbinom_sample)
hist(nbinom_sample)
write.csv(nbinom_sample, file = "C:/Users/Delaram/Desktop/coverage_sample.csv")

sim.read.count = sum(nbinom_sample)

