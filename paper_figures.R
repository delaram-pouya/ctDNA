## plotting the figures needed for paper


source('error/Codes/4_findError_functions.R')
Initialize()

############# mixture guassian cfDNA error 
####### large peak 

large_index = findLargePeakIndex(mixEM_list)

mu_large_peak <- sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$mu[ large_index[i]  ]  )
sigma_large_peak <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$sigma[ large_index[i]  ]  )
pmix_large_peak <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$lambda[ large_index[i]  ]  )


###### small peak

small_index = findSmallPeakIndex(mixEM_list)

mu_small_peak <- sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$mu[ small_index[i]  ]  )
sigma_small_peak <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$sigma[ small_index[i]  ]  )
pmix_small_peak <-  sapply(1:length(mixEM_list), function(i) mixEM_list[[i]]$lambda[ small_index[i]  ]  )


names(mu_large_peak)= Names; names(sigma_large_peak) = Names; names(pmix_large_peak) = Names
names(mu_small_peak)= Names; names(sigma_small_peak) = Names; names(pmix_small_peak) = Names

###########

real_data = data.frame(ListOfSplitednonZeroErrorsMergedLog[['cfDNA']])
colnames(real_data) <- 'empirical_error'

n = nrow(real_data)
data_large_peak = data.frame(rnorm(round(n*pmix_large_peak[['cfDNA']]), mu_large_peak[['cfDNA']], sigma_large_peak[['cfDNA']] ))
data_small_peak = data.frame(rnorm(round(n*pmix_small_peak[['cfDNA']]), mu_small_peak[['cfDNA']], sigma_small_peak[['cfDNA']] ))

colnames(data_large_peak) <- 'estimated_error'
colnames(data_small_peak) <- 'estimated_error'


ggplot(real_data, aes(x=empirical_error))+geom_density(aes(y = ..count..),color='white',fill='blue', alpha=0.2)+
  geom_density(data=data_large_peak,aes(x=estimated_error, y = ..count..), color='dark blue',linetype='dashed')+ stat_density(geom="line")+
  geom_density(data=data_small_peak,aes(x=estimated_error, y = ..count..), color='red',linetype='dashed')+ stat_density(geom="line")+
  theme_bw()



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

P_error <- ggplot(real_data, aes(x=empirical_error))+geom_density(aes(y = ..count.., fill='empirical'),color='black', alpha=0.4)+
  geom_density(data=data_large_peak,aes(x=estimated_error, y = ..count.., colour='large-peak'), linetype='solid',size=0.95)+ stat_density(geom="line")+
  geom_density(data=data_small_peak,aes(x=estimated_error, y = ..count.., colour='small-peak'),linetype='solid', size=0.95)+ stat_density(geom="line")+
  theme_classic()+
  scale_colour_manual(name='Estimated Error' ,values = c('#E69F00', '#0072B2'))+
  scale_fill_manual(name='Empirical Error', labels='',values = '#999999')+
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 2, override.aes = list(alpha = 0.5,colour='black',linetype='solid',size=0.7)))+
  theme(legend.position = c(0.8, 0.5))+xlab('cfDNA Error (log scaled)')


pdf('P_error_mix_cfDNA.pdf', height=7, width=8)
P_error
dev.off()





#############################################

## Tumor origin distribution

pileupVAF <- readRDS('simulator/Tumor_origin_read_probability/pileUpVAFall.rds')
vcf <- readRDS('simulator/Tumor_origin_read_probability/allVCFs.rds')


## merging data
merged_VAFs <- c(pileupVAF[['p1_cfDNA1']]$VAF,
                 pileupVAF[['p1_cfDNA2']]$VAF, 
                 pileupVAF[['p2_cfDNA1']]$VAF,
                 pileupVAF[['p2_cfDNA2']]$VAF)


merged_VAFs <- data.frame(merged_VAFs[merged_VAFs>0])
colnames(merged_VAFs) <- 'VAF'


#### generate samples from all distributions
fitParameters <- readRDS('simulator/Tumor_origin_read_probability/fitParameters_lev2_corrected.rds')
SAMPLE_SIZE = nrow(merged_VAFs)

beta <- data.frame(rbeta(SAMPLE_SIZE, fitParameters$beta$estimate[1],fitParameters$beta$estimate[2]))
normal <- data.frame(rnorm(SAMPLE_SIZE, fitParameters$norm$estimate[1],fitParameters$norm$estimate[2]))
gamma <- data.frame(rgamma(SAMPLE_SIZE, fitParameters$gamma$estimate[1],fitParameters$gamma$estimate[2]))
lnorm <- data.frame(rlnorm(SAMPLE_SIZE, fitParameters$lnorm$estimate[1],fitParameters$lnorm$estimate[2]))
weibull <- data.frame(rweibull(SAMPLE_SIZE, fitParameters$weibull$estimate[1],fitParameters$weibull$estimate[2]))

colnames(beta) <- 'beta'
colnames(normal) <- 'normal'
colnames(gamma) <- 'gamma'
colnames(lnorm) <- 'lnorm'
colnames(weibull) <- 'weibull'



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf('P_tissueProb_cfDNA.pdf')
ggplot(merged_VAFs, aes(x=VAF))+geom_density(aes(fill='empirical'),color='yellow',size=0.1, alpha=0.3)+
  
  geom_density(data=normal,aes(x=normal, colour='normal'), linetype='solid',size=1)+ stat_density(geom="line")+
  geom_density(data=gamma,aes(x=gamma, colour='gamma'),linetype='solid', size=1)+ stat_density(geom="line")+
  geom_density(data=lnorm,aes(x=lnorm, colour='lnorm'),linetype='solid', size=1)+ stat_density(geom="line")+
  geom_density(data=beta,aes(x=beta, colour='beta'),linetype='solid', size=1)+ stat_density(geom="line")+
  #geom_density(data=weibull,aes(x=weibull, colour='weibull'),linetype='solid', size=1)+ stat_density(geom="line")+
  
  theme_classic()+
  scale_colour_manual(name='Estimated P_tis' ,values = c( "#CC79A7", "#009E73" ,"#0072B2", "#D55E00"))+
  scale_fill_manual(name='Empirical P_tis', labels='',values = 'gold1')+
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 2, override.aes = list(alpha = 0.5,colour='#E69F00',linetype='solid',size=0.7)))+
  theme(legend.position = c(0.77, 0.5))+xlab('cfDNA P_tis')+xlim(c(0,1))
dev.off()


################################



depth <- readRDS('simulator/winow_depth/all_chr_converage.rds')
poisson_parameters <- fitdist(data = depth, distr = 'pois', discrete = T)
nbinomial_parameters <- fitdist(data = depth, distr = 'nbinom', discrete = T)

n = length(depth)
pois <- data.frame(rpois(n, poisson_parameters$estimate[1]))
nbinom <- data.frame(rnbinom(n, size=nbinomial_parameters$estimate[1], mu=nbinomial_parameters$estimate[2]))
colnames(pois) <- 'pois'
colnames(nbinom) <- 'nbinom'

summary(depth)
summary(pois)
summary(nbinom)

ggplot(data.frame(depth), aes(x=depth))+geom_histogram(color='black', fill='pink', alpha=0.8,bins = 20)+
  xlim(c(1,200))+theme_classic()

ggplot(data.frame(depth), aes(x=depth)) +
  geom_histogram(color='black', fill='pink', alpha=0.8,breaks=c(seq(0, 200, by=20), max(depth)), position = "identity") +
  coord_cartesian(xlim=c(0,220))


ggplot(data.frame(depth), aes(x=depth))+
  geom_histogram(color='darkblue', fill='lightblue', alpha=0.6,breaks=c(seq(0, 200, by=20)))+
  scale_x_continuous(limits=c(0, 200), breaks=c(seq(0, 200, by=20)), labels=c(seq(0,190, by=20), "200<"))+
  theme_classic()+xlab('window-depth')

pdf('P_coverage_cfDNA.pdf')
ggplot(data.frame(depth), aes(x=depth))+
  geom_density(aes(fill='empirical'),color='gray40', alpha=0.7)+
  geom_density(data=nbinom, aes(nbinom,fill='negative binomial'), color='black', alpha=0.6,size=0.4)+
  theme_classic()+xlab('window-coverage')+
  scale_fill_manual(name='',values = c('lightslategray', 'red2'))+
  guides(fill = guide_legend(override.aes = list(alpha = 0.5)))+
  theme(legend.position = c(0.77, 0.5))+
  scale_x_continuous(limits=c(0, 100), breaks=c(seq(0, 100, by=50)), labels=c(seq(0,50, by=50), "100<"))
dev.off()





