## In this script, we will score the genes  

source('annotate/Codes/1_Functions.R')
Initialize()

#### importing datasets needed for scoring and annotation
CGC <- read.csv('annotate/Sources/CancerResources/Cancer_Gene_Census.csv')
DGIdb <- read.csv('annotate/Sources/CancerResources/DGIdb_interactions.csv')
PanCan_Sig <- read.csv('annotate/Sources/CancerResources/PanCan_cancertype.csv')
PanCan_Driver <- read.csv('annotate/Sources/CancerResources/PanCan_drivers.csv')
oncoKB <- read.delim('annotate/Sources/CancerResources/oncoKB_allCuratedGenes.txt')
listOfDataSets <- list(CGC, DGIdb, PanCan_Driver, oncoKB)
names(listOfDataSets) <- c('CGC', 'DGIdb', 'PanCan_Driver', 'oncoKB')





qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

## takes the annotated vcf file as the input 
FinalAnnotatedVcf <- read.csv('Annaotated_VCF.csv')


## make gene table
geneSet <- FinalAnnotatedVcf$hgnc_symbol
geneSet[is.na(geneSet)|geneSet=='']='Not_Annotated'
geneSetDf <- data.frame(table(geneSet))
rownames(geneSetDf) <- geneSetDf$geneSet


listOfDataSetGenes <- lapply(listOfDataSets, function(x) toupper(x$Gene))
isGeneSetIncluded <- data.frame(lapply(listOfDataSetGenes, function(x) ifelse( geneSetDf$geneSet %in% x, 1, 0)  ))

ScoreTable <- cbind(geneSetDf, isGeneSetIncluded)
ScoreTable <- ScoreTable[rownames(ScoreTable) != 'Not_Annotated',]


### Visualize
pdf('HighScoreGenes_barplots.pdf')
if( isEmpty(geneSetDf) ){
    grid.text("No annotated gene in the vcf file ",x = 0.5, y=0.5 ,gp=gpar(fontsize=12),just = "center")
    grid.newpage()}


ScoreTable$geneSet <- factor(ScoreTable$geneSet, levels = ScoreTable$geneSet[order(ScoreTable$Freq)])
ggplot(subset(ScoreTable,Freq>quantile(Freq, 0.75)), aes(x=geneSet, y=Freq, fill=Freq))+ggtitle('Highly Mutated Genes')+
    geom_bar(stat="identity",color="black",alpha=0.85)+coord_flip()+theme_bw()+
    scale_fill_gradient2(low = "white", high = "blue")+xlab('Gene Set')+ylab('Mutation Frequency')


ScoreTable$normalized_Freq <- getNormalizeFreq(ScoreTable$Freq)
ScoreTable$Score <- rowSums(subset(ScoreTable, select=-c(geneSet, Freq)))
ScoreTable <- ScoreTable[order(ScoreTable$Score, decreasing = T),]
head(ScoreTable)
write.csv(ScoreTable, 'ScoreTable.csv')



ScoreTable$geneSet <- factor(ScoreTable$geneSet, levels = ScoreTable$geneSet[order(ScoreTable$Score)])
data <- subset(ScoreTable,Score>quantile(Score, 0.75))
ggplot(data , aes(x=geneSet, y=Score,fill=geneSet))+ggtitle('High-score Genes')+
    geom_bar(stat="identity",alpha=0.9,color='black')+coord_flip()+theme_bw()+ theme(legend.position = "none") +
  scale_fill_manual(values = sample(col_vector, nrow(data)))+xlab('Gene Set')
dev.off()
