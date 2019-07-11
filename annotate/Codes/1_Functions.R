# function to check to see if packages are installed. 
#  Install them if they are not, then load them into the R session.

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}




## loading required packages    

Initialize <- function(){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    options(stringsAsFactors = FALSE)
    listOfPackages <- c('Biobase',  'limma',  'pheatmap', 'gplots', 'reshape2', 
                        'plyr', 'stringr', 'ggplot2',  'parallel', "AnnotationDbi", "grid",
                        "gridExtra", "magrittr","dplyr","pillar" , "plyr" , "Rcpp", 
                        'VariantAnnotation', 'tidyr', 'org.Hs.eg.db', 'reactome.db', 'clusterProfiler',
                        'ReactomePA', 'vcfR', 'data.table', 'biomaRt')
    ipak(listOfPackages)
}



## importing datasets 

ImportDataSets <- function(){
    
    SourceFiles <- list.files('annotate/Sources/CancerResources', full.names = T, pattern = '*.csv')
    SourceNames <- gsub('.csv','',gsub('annotate/Sources/CancerResources/', '', SourceFiles))
    SourceDataSets <- lapply(SourceFiles, read.csv, header= T)
    names(SourceDataSets) <- SourceNames
    
    # literature <- data.frame(c("ESR1","ERBB2","PIK3CA","AKT1","GATA3","HER2", "TP53","BRCA1", "BRCA2"))
    # colnames(literature) <- c("Gene")
    # SourceDataSets <- append(list(literature), SourceDataSets, 1)
    # names(SourceDataSets)[1] <- 'literature'
    
    names(SourceDataSets) <- c( 'oncoKB_actioable', 'CGC',
                               'DGIdb', 'PanCan', 'DrugBank', 'oncoKB', 'PanCanDriver', 'PanCanSig')
    return(SourceDataSets)
} 



isEmpty <- function(table){ ifelse(nrow(table) == 0, TRUE, FALSE)}


cleanTranscriptId <- function(selectedCol){
    selectedCol_split <- str_split(selectedCol$V1, ' ')
    transcriptID <- lapply(selectedCol_split, function(x) x[length(x)])
    gsub('[\"|;]', '', transcriptID, perl=T)
}




getCleanedAnnotationFile <- function(hg38Annotations){
    hg38Annotations = fread('annotate/Sources/hg38GeneAnnotations.txt', header = F)
    hg38Annotations$transcriptId <- cleanTranscriptId(data.table(hg38Annotations$V9))
    hg38Annotations <- subset(hg38Annotations, select=-c(V2,V9,V6,V7,V8))
    colnames(hg38Annotations)[-ncol(hg38Annotations)] <- c('chr','region','start','end')
    return(hg38Annotations)
}



getNormalizeFreq <- function(x){
  ifelse(x>quantile(x, 0.75), 3, ifelse(x>median(x), 2, ifelse(x>quantile(x, 0.25), 1, 0))  ) 
}


