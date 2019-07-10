## In this script, we will 

source('Codes/1_Functions.R')
Initialize()

## takes the input directory and vcf file name as an input 
input_DIR = '~/Desktop/'
vcf_name = 'p2_cfDNA1_snv_lev3.vcf'

## makes a directory to write the outputs in it
outout_DIR = paste0(input_DIR, gsub('.vcf','',vcf_name),'/')
output_prefix = paste0(outout_DIR, gsub('.vcf','',vcf_name))

FinalAnnotatedVcf <- read.csv(paste0(output_prefix,'_annaotatedVCF.csv'))


## make gene table
geneSet <- FinalAnnotatedVcf$hgnc_symbol
geneSet[is.na(geneSet)|geneSet=='']='NotAnnot'
geneSetDf <- data.frame(table(geneSet))
rownames(geneSetDf) <- geneSetDf$geneSet



listOfDataSets <- ImportDataSets() 
listOfDataSetGenes <- lapply(listOfDataSets, function(x) toupper(x$Gene))
isGeneSetIncluded <- data.frame(lapply(listOfDataSetGenes, function(x) ifelse( geneSetDf$geneSet %in% x, 1, 0)  ))


ScoreTable <- cbind(geneSetDf, isGeneSetIncluded)
ScoreTable$Score <- rowSums(subset(ScoreTable, select=-geneSet))
ScoreTable <- ScoreTable[order(ScoreTable$Score, decreasing = T),]
head(ScoreTable)
write.csv(ScoreTable, paste0(output_prefix, '_ScoreTable.csv'))


### Visualize
pdf(paste0(output_prefix,'_highScoreGenes_barplots.pdf'))
if( isEmpty(geneSetDf) ){
    grid.text("No annotated gene in the vcf file ",x = 0.5, y=0.5 ,gp=gpar(fontsize=12),just = "center")
    grid.newpage()}
ggplot(subset(ScoreTable,Score>mean(Score)), aes(x=geneSet, y=Freq))+ggtitle('Highly mutated genes')+
    geom_bar(stat="identity",color="black", fill='#E69F00',alpha=0.6)+coord_flip()+theme_bw()
ggplot(subset(ScoreTable,Score>mean(Score)), aes(x=geneSet, y=Score))+ggtitle('High-score genes')+
    geom_bar(stat="identity",color="black", fill='cyan3',alpha=0.6)+coord_flip()+theme_bw()
dev.off()

       
