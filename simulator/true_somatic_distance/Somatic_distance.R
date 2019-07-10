## in this script we'll try to find the 
##  distribution of sequential somatic mutations' distance from each other  

options(stringsAsFactors = F)
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


## import data
trueSomatics <- read.csv("Simulator/true_somatic_distance/p1_true_somatics.txt", header = F, stringsAsFactors = F)
trueSomatics <- as.data.frame(str_split_fixed(trueSomatics$V1, "_", 2))
colnames(trueSomatics) <- c("chr", "position")
head(trueSomatics)


chromosomes <- paste0('chr',c(1:22,'X', 'Y'))


listOfDistances <- lapply(chromosomes, function(CHROM){
  ## split based on chr
  chrom_positions <- subset(trueSomatics, chr==CHROM)
  chrom_positions$position <- as.numeric(chrom_positions$position)
  ## sort the positions
  chrom_positions <- chrom_positions[order(chrom_positions$position),]
  ## find the distance between positions
  distace <- unlist(sapply(1:length(chrom_positions$position), function(i) {
    if( i%%2==0) chrom_positions$position[i]-chrom_positions$position[i-1]}) )
  ## if number of positions is odd, add the last distance
  if(nrow(chrom_positions)%%2==1) distace = c(distace, 
                                              chrom_positions$position[nrow(chrom_positions)] - chrom_positions$position[nrow(chrom_positions)-1])
  if(length(distace)) return(distace)
})

distances <- unlist(listOfDistances)
ggplot(data.frame(distances), aes(x=distances))+geom_density()+theme_bw()+ggtitle('true somatic distance distribution')

sum(distances<300)/length(distances) ## only 0.1 of distances are bellow 300 



    