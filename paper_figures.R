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
  geom_density(data=data_large_peak,aes(x=estimated_error, y = ..count.., colour='large-peak'), linetype='dashed',size=0.95)+ stat_density(geom="line")+
  geom_density(data=data_small_peak,aes(x=estimated_error, y = ..count.., colour='small-peak'),linetype='dashed', size=0.95)+ stat_density(geom="line")+
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

lapply(pileupVAF, head)
lapply(vcf, head)

## merging data
merged_VAFs <- c(pileupVAF[['p1_cfDNA1']]$VAF,
                 pileupVAF[['p1_cfDNA2']]$VAF, 
                 pileupVAF[['p2_cfDNA1']]$VAF,
                 pileupVAF[['p2_cfDNA2']]$VAF)


merged_VAFs <- data.frame(merged_VAFs[merged_VAFs>0])
colnames(merged_VAFs) <- 'VAF'

ggplot(data=merged_VAFs, aes(VAF)) + 
  geom_density(col="blue", fill="#56B4E9", alpha = .2) + 
  labs(title="tumor-origin reads distribution") +
  labs(x=" %confirming alterations", y="frequency") +theme_bw()


fitParameters <- readRDS('simulator/Tumor_origin_read_probability/fitParameters_lev2_corrected.rds')



# create some data to work with
x = rnorm(1000);

# overlay histogram, empirical density and normal density
p0 = qplot(x, geom = 'blank') +   
  geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Density', values = c('red', 'blue')) + 
  theme(legend.position = c(0.85, 0.85))
p0

