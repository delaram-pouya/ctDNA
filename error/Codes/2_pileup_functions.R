## This script includes the functions which are used in the pileup.R file

Initialize <- function(){
  options(stringsAsFactors = F)
  library(plyr)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(ggplot2)
  library(Rsamtools)
  library(vcfR)
  library(gridExtra)
}


MakepileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tablecounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  return(res)}


Add.attrib <- function(df){
  df$depth <- df$A + df$T + df$C + df$G
  tmp <- t(apply(subset(df,select=c('A','C','G','T')), 1, sort))
  
  df$F1 <- tmp[,4]/df$depth # Maximum frequency allele
  df$F2 <- tmp[,3]/df$depth # Second maximum frequency allele
  df$ALT.freq <- (rowSums(tmp[,-4]))/df$depth # (rowSum - max1)/depth
  head(df)
  df.n <- subset(df, select=c('A','T','C','G'))
  m1.index <- apply(df.n,1,which.max)
  df$m1 <- colnames(df.n)[apply(df.n,1,which.max)]
  
  sapply(1:nrow(df.n), function(i) {df.n[i,m1.index[i]] <<-(-1)}) # faster than sorting
  df$m2 <- colnames(df.n[,1:4])[apply(df.n[,1:4],1,which.max)]
  df[(is.na(df$ALT.freq))|(df$ALT.freq==0),]$m2<-NA
  
  return(df)}



## keeping only the position which are present in all files of a patient
find.common <- function(p){
  id <-lapply(p, function(x) paste0(x$seqnames,sep='_',x$pos) )
  p.common <- Reduce(intersect,id)
  p <- sapply(1:length(p), function(i) p[[i]][id[[i]]%in%p.common,] , simplify = F)
  return(p)}




### calculating VAF
find.VAF <- function(i,df){
  x=df[i,]
  if(is.na(x$ALT) | x$depth==0 | nchar(x$ALT)>1) return(0)
  else if(!is.na(x$ALT)){base=x$ALT;return(x[[base]]/x$depth)}
}




