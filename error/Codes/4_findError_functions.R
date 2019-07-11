## this scipt includes the functions needed in the 
## findError.R file


Initialize <- function(){
  options(stringsAsFactors = F)
  library(dplyr)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(ggplot2)
  library(Rsamtools)
  library(vcfR)
  library(gridExtra)
  library(reshape2)
  library(fitdistrplus) 
  require(MASS) 
  library(mixtools)
}



preprop <- function(df){
  df <- subset(df, seqnames %in% paste0('chr',c(1:22,c('X','Y'))) )
  df <- subset(df, select=c('seqnames','pos','nucleotide','count'))
  df <- subset(df, nucleotide %in% c('A','T','C','G'))
  return(df)
}



MakepileupFreq.eff <- function(df){
  df = df[, sum(count), by = list(seqnames, pos, nucleotide)]
  setnames(df, 'V1', 'sumCount')
  df = dcast.data.table(df, seqnames + pos ~ nucleotide, sum, value.var = 'sumCount')
  return(df)
}



Add.depth <- function(df){
  df$depth <- df$A + df$T + df$C + df$G
  return( subset(df, depth>20) ) #depth-filter=20
}



## homo-filter 
Add.max.blood <- function(df){ 
  df <- df[, mVal := do.call(pmax, .SD), .SDcols = 3:6]
  #df <- df[,mVal:=max(A, T, C, G), by=1:nrow(df)] #alternative 
  df$F1 <- df$mVal /df$depth
  df$ALT.freq <- 1 - df$F1
  df <-  df[,m1 :=  names(.SD)[max.col(.SD)], .SDcols = 3:6]
  return(df[ df$F1>0.8 ,!c('mVal')]) ## homo-thershold = 0.8
}




## No homo-filter 
Add.max <- function(df){
  df <- df[, mVal := do.call(pmax, .SD), .SDcols = 3:6]
  #df <- df[,mVal:=max(A, T, C, G), by=1:nrow(df)] #alternative 
  df$F1 <- df$mVal /df$depth
  df$ALT.freq <- 1 - df$F1
  df <-  df[,m1 :=  names(.SD)[max.col(.SD)], .SDcols = 3:6]
  return(df[,!c('mVal')]) 
}




matUpdate <- function(mat,x){
  sapply(1:nrow(x),function(i){ 
    row <- x[i,]
    max = row$m1
    alt <- colnames(row[,1:4])[colnames(row[,1:4])!=max] 
    altFreq <- (as.numeric(row[,1:4])[colnames(row[,1:4])!=max] )/row$depth #count of bases other than F1
    for(t in 1:length(alt)){
      mut <- alt[t] 
      mat[[max,mut]] <<- c( mat[[max,mut]],  altFreq[t] ) }
  })
  return(mat)
} 




makeErrorTable <- function(dataMat){
  A <-  dataMat[1,];T <- dataMat[2,]; C <- dataMat[3,]; G <- dataMat[4,]
  
  Amut <- data.frame(val=c(A$T[[1]],A$C[[1]],A$G[[1]]),
                     var=c( rep('A_T',length(A$T[[1]])), 
                            rep('A_C',length(A$C[[1]])), 
                            rep('A_G',length(A$G[[1]])) ) )
  
  Tmut <- data.frame(val=c(T$A[[1]],T$C[[1]],T$G[[1]]),
                     var=c( rep('T_A',length(T$A[[1]])), 
                            rep('T_C',length(T$C[[1]])), 
                            rep('T_G',length(T$G[[1]])) ) )
  
  Cmut <- data.frame(val=c(C$A[[1]],C$T[[1]],C$G[[1]]),
                     var=c( rep('C_A',length(C$A[[1]])), 
                            rep('C_T',length(C$T[[1]])), 
                            rep('C_G',length(C$G[[1]])) ) )
  
  Gmut <- data.frame(val=c(G$A[[1]],G$T[[1]],G$C[[1]]),
                     var=c( rep('G_A',length(G$A[[1]])), 
                            rep('G_T',length(G$T[[1]])), 
                            rep('G_C',length(G$C[[1]])) ) )
  
  All.mut <- do.call(rbind, list(Amut,Tmut,Cmut,Gmut))
  return(All.mut)
}


GetBaseWiseError <- function(dataMat){
  A <-  dataMat[1,];T <- dataMat[2,]; C <- dataMat[3,]; G <- dataMat[4,]
  mutationNames = c('A_T','A_C','A_G',  'T_A','T_C','T_G',  'C_A','C_T','C_G',  'G_A','G_T','G_C')
  mutationList   <- list(  A$T[[1]], A$C[[1]], A$G[[1]],  
                      T$A[[1]], T$C[[1]], T$G[[1]],  
                      C$A[[1]], C$T[[1]], C$G[[1]],
                      G$A[[1]], G$T[[1]], G$C[[1]])
  names(mutationList) = mutationNames
  return(mutationList)
}


makePlot <- function(All.mut, FileName){
  
  p1 <- ggplot(All.mut,aes(val+1e-7, color=var))+geom_density()+scale_x_log10()+theme_bw()+
    scale_color_brewer(palette = 'Paired')+xlab('alteration rate(log)')+ 
    ggtitle(paste0(FileName,' nucleotide-wise error Rate(homozygous)'))
  
  p2 <- ggplot(All.mut,aes(y=val+1e-7, x=var))+geom_violin(aes(fill=var))+
    scale_y_log10()+theme_bw()+scale_fill_brewer(palette = 'Paired')+ylab('alteration rate(log)')
  
  #nonZero <- lapply(list(Amut,Tmut,Cmut,Gmut), subset, val>0)
  #all.mut.nonZ <- do.call(rbind, nonZero)
  all.mut.nonZ <- subset(All.mut, val>0)
  
  p3 <- ggplot(all.mut.nonZ,aes(val, color=var))+geom_density()+scale_x_log10()+
    theme_bw()+scale_color_brewer(palette = 'Paired')+
    ylab('alteration rate(log)')+ 
    ggtitle(paste0(FileName,' nucleotide-wise error Rate(homozygous & non-Zero)'))
  
  
  p4 <- ggplot(all.mut.nonZ,aes(x=var,y=val))+geom_violin(aes(fill=var))+
    scale_y_log10()+theme_bw()+scale_fill_brewer(palette='Paired')+ylab('alteration rate(log)')
  
  pdf(paste0(FileName,'_base_mutations.pdf'),height=15,width=18)
  gridExtra::grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
  dev.off()
  
}



FindZeroRatio <- function(x){
  SampleRatios <- data.frame(table(x$val==0))
  rownames(SampleRatios) <- NULL
  colnames(SampleRatios) <- c('variable','Freq')
  SampleRatios$variable <- c('non-Zero','Zero')
  SampleRatios$Freq <- round(SampleRatios$Freq /sum(SampleRatios$Freq ),2)
  return(SampleRatios)
}


AddAttributes <- function(x){
  x <- mclapply(x, preprop, mc.cores = detectCores()-2)
  x <- mclapply(x, MakepileupFreq.eff, mc.cores = detectCores()-2)
  x <- mclapply(x, Add.depth , mc.cores = detectCores()-2)
  x <- mclapply(x, Add.max , mc.cores = detectCores()-2) 
  return(x)}


AddAttributes_singleDF <- function(x){
  x <- preprop(x)
  x <- MakepileupFreq.eff(x)
  x <- Add.depth(x)
  x <- Add.max(x)
  return(x)}



findLargePeakIndex <- function(mixEM_list){
  minMu <- lapply(mixEM_list, function(x) min(x$mu[1], x$mu[2]))
  sapply(1:length(minMu), function(i)ifelse(mixEM_list[[i]]$mu[1] ==minMu[[i]], 1 , 2) )
}


findSmallPeakIndex <- function(mixEM_list){
  maxMu <- lapply(mixEM_list, function(x) max(x$mu[1], x$mu[2]))
  sapply(1:length(maxMu), function(i)ifelse(mixEM_list[[i]]$mu[1] ==maxMu[[i]], 1 , 2) )
}


