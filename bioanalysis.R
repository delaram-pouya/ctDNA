#!/usr/bin/env Rscript

## This script, gets a vcf file as an input and annotates it gencode-v29(hg38)
## reference genome, then we'll integrate the input data with some datasets to
## find the most important genes
## NOTE: MAKE SURE THE COORDINATIONS OF VCF FILE IS BASED ON HG38 REFERENCE GENOME

## input: vcf file
## output: 
# 1. annotated csv file with gene symboles
# 2. score table csv file
# 3. pdf file (analysis results + pathway enrichment)

# run as:
# Rscript bioanalysis.R annotate/Data/p2_cfDNA1_snv_lev3.vcf annotate/Sources/hg38GeneAnnotations.txt annotate/Sources/CancerResources/Cancer_Gene_Census.csv annotate/Sources/CancerResources/DGIdb_interactions.csv annotate/Sources/CancerResources/PanCan_cancertype.csv annotate/Sources/CancerResources/PanCan_drivers.csv annotate/Sources/CancerResources/oncoKB_allCuratedGenes.txt

# args[1]  = vcf file
# args[2]  = gene-code annotations (hg38GeneAnnotations.txt)
# args[3]  = Cancer_Gene_Census.csv
# args[4]  = DGIdb_interactions.csv
# args[5]  = PanCan_cancertype.csv
# args[5]  = PanCan_drivers.csv
# args[5]  = oncoKB_allCuratedGenes.txt


########################## Functions
args = commandArgs(trailingOnly=TRUE)


ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    BiocManager::install(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


Initialize <- function(){
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  options(stringsAsFactors = FALSE)
  listOfPackages <- c('Biobase', 'gplots', 'reshape2', 'RColorBrewer',
                      'plyr', 'stringr', 'ggplot2',  'parallel', "AnnotationDbi", "grid",
                      "gridExtra", "magrittr","dplyr","pillar" , "plyr" , "Rcpp", 
                      'VariantAnnotation', 'tidyr', 'org.Hs.eg.db', 'reactome.db', 'clusterProfiler',
                      'ReactomePA', 'vcfR', 'data.table', 'biomaRt')
  ipak(listOfPackages)
}


isEmpty <- function(table){ ifelse(nrow(table) == 0, TRUE, FALSE)}


cleanTranscriptId <- function(selectedCol){
  selectedCol_split <- str_split(selectedCol$V1, ' ')
  transcriptID <- lapply(selectedCol_split, function(x) x[length(x)])
  gsub('[\"|;]', '', transcriptID, perl=T)
}


getCleanedAnnotationFile <- function(hg38Annotations){
  hg38Annotations = fread(args[2], header = F)
  hg38Annotations$transcriptId <- cleanTranscriptId(data.table(hg38Annotations$V9))
  hg38Annotations <- subset(hg38Annotations, select=-c(V2,V9,V6,V7,V8))
  colnames(hg38Annotations)[-ncol(hg38Annotations)] <- c('chr','region','start','end')
  return(hg38Annotations)
}


getNormalizeFreq <- function(x){
  ifelse(x>quantile(x, 0.75), 3, ifelse(x>median(x), 2, ifelse(x>quantile(x, 0.25), 1, 0))  ) 
}

Initialize()

## define color palette for visualization 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


########################## Importing inputs 

## takes the vcf file as an input 
VCF_INPUT = args[1]

## import gencode-v29(hg38) annotation file
AnnotationHg38 <- getCleanedAnnotationFile()

#### importing datasets needed for scoring and annotation
CGC <- read.csv(args[3])
DGIdb <- read.csv(args[4])
PanCan_Sig <- read.csv(args[5])
PanCan_Driver <- read.csv(args[6])
oncoKB <- read.delim(args[7])

listOfDataSets <- list(CGC, DGIdb, PanCan_Driver, oncoKB)
names(listOfDataSets) <- c('CGC', 'DGIdb', 'PanCan_Driver', 'oncoKB')



## import vcf file
InputVcf <- read.vcfR( VCF_INPUT, verbose = FALSE )
InputVcf <- vcfR2tidy(InputVcf, info_only = T, single_frame = FALSE, toss_INFO_column = TRUE)
InputVcf_fix <- (InputVcf$fix)
InputVcf_fix <- InputVcf_fix[rowSums(is.na(InputVcf_fix)) != ncol(InputVcf_fix), ]  ## removing rows with all NA values


### annotating vcf with gencode-v29 hg38 
query = GRanges(seqnames = InputVcf_fix$CHROM, IRanges(start = InputVcf_fix$POS, end=InputVcf_fix$POS, width=1))
subject = GRanges(AnnotationHg38$chr, IRanges(start = AnnotationHg38$start, end=AnnotationHg38$end))
hits <- findOverlaps(query, subject, type='within')
VcfAnnaotatedByEnsembl <- cbind(data.frame(InputVcf_fix[queryHits(hits),]), data.frame(AnnotationHg38[subjectHits(hits),]))


## removing version from ensembl transcript IDs
VcfAnnaotatedByEnsembl$transcriptId <- gsub('\\..*' ,'',VcfAnnaotatedByEnsembl$transcriptId)


## converting ensmble IDs to hgnc_symbol and entrezgene
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
GeneSymbols <- getBM(filters= "ensembl_transcript_id", 
                     attributes= c("ensembl_transcript_id", "hgnc_symbol","description"),
                     values= VcfAnnaotatedByEnsembl$transcriptId,
                     mart= mart)


FinalAnnotatedVcf <- merge(VcfAnnaotatedByEnsembl, GeneSymbols, by.x='transcriptId', by.y= 'ensembl_transcript_id',all.x=T)
head(FinalAnnotatedVcf)

entrezgene <- bitr(FinalAnnotatedVcf$hgnc_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
colnames(entrezgene) <- c('SYMBOL', 'entrezgene')
FinalAnnotatedVcf <- merge(FinalAnnotatedVcf,entrezgene, by.x='hgnc_symbol', by.y='SYMBOL',all.x=TRUE)
head(FinalAnnotatedVcf)

write.csv(FinalAnnotatedVcf, file ='Annaotated_VCF.csv', row.names = F)

########################## Score genes

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

ScoreTable$geneSet <- factor(ScoreTable$geneSet, levels = ScoreTable$geneSet[order(ScoreTable$Freq)])
p_mutFreq <- ggplot(subset(ScoreTable,Freq>quantile(Freq, 0.75)), aes(x=geneSet, y=Freq, fill=Freq))+ggtitle('Highly Mutated Genes')+
  geom_bar(stat="identity",color="black",alpha=0.85)+coord_flip()+theme_bw()+
  scale_fill_gradient2(low = "white", high = "blue")+xlab('Gene Set')+ylab('Mutation Frequency')


ScoreTable$normalized_Freq <- getNormalizeFreq(ScoreTable$Freq)
ScoreTable$Score <- rowSums(subset(ScoreTable, select=-c(geneSet, Freq)))
ScoreTable <- ScoreTable[order(ScoreTable$Score, decreasing = T),]
head(ScoreTable)
write.csv(ScoreTable, 'ScoreTable.csv')



ScoreTable$geneSet <- factor(ScoreTable$geneSet, levels = ScoreTable$geneSet[order(ScoreTable$Score)])
data <- subset(ScoreTable,Score>quantile(Score, 0.75))
p_geneScore <- ggplot(data , aes(x=geneSet, y=Score))+ggtitle('High-score Genes')+xlab('Gene Set')+
  geom_bar(stat="identity",alpha=0.75,color='black',fill='pink2')+coord_flip()+theme_bw()+ theme(legend.position = "none") 



########################## Gene Annotation table

lapply(listOfDataSets, head)
listOfDataSetGenes <- lapply(listOfDataSets, function(x) toupper(x$Gene))
colnames(geneSetDf)[1] <- c('Gene')

## preprocess CGC
CGC <- subset(CGC, Gene %in% geneSet)
CGC <- subset(CGC, select=c(Gene, Name, Tier, Chr.Band, Tumour.Types.Germline. ,Tumour.Types.Somatic., 
                            Germline,Somatic, Role.in.Cancer, Mutation.Types ,Translocation.Partner))
colnames(CGC)[2:ncol(CGC)] <- paste0(colnames(CGC)[2:ncol(CGC)] , '(CGC)')
head(CGC)


## PanCan_Sig
PanCan_Sig <- subset(PanCan_Sig, Gene %in% geneSet )
PanCan_Sig <- subset(PanCan_Sig, select=c(Gene, Pathway) )
colnames(PanCan_Sig)[2:ncol(PanCan_Sig)] <- paste0(colnames(PanCan_Sig)[2:ncol(PanCan_Sig)] , '(PanCan)') 
head(PanCan_Sig)    


## PanCan drivers
colnames(PanCan_Driver)[3] <- 'tsg_onco_20.20pred'
PanCan_Driver <- subset(PanCan_Driver, Gene %in% geneSet)
PanCan_Driver <- subset(PanCan_Driver, select=c(Gene, Cancer, tsg_onco_20.20pred,Pancan.Frequency))
colnames(PanCan_Driver)[2:ncol(PanCan_Driver)] <- paste0(colnames(PanCan_Driver)[2:ncol(PanCan_Driver)] , '(PanCan)')
head(PanCan_Driver)


## DGIdb 
DGIdb <- subset(DGIdb, Gene %in% geneSet)
DGIdb <- subset(DGIdb, select=c(Gene, drug_name, drug_chembl_id))
colnames(DGIdb)[2:ncol(DGIdb)] <- paste0(colnames(DGIdb)[2:ncol(DGIdb)] , '(DGIdb)')
head(DGIdb)


merged_table = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Gene", all.x = TRUE),list(geneSetDf,PanCan_Driver,PanCan_Sig,CGC,DGIdb))
merged_table <- merged_table[order(merged_table$Freq, decreasing = T),]
merged_table <- merged_table[merged_table$Gene != 'Not_Annotated',]
rownames(merged_table) <- NULL
head(merged_table)
write.csv(merged_table, 'GeneInfoTable.csv')


########################## Pathway Enrichment

entrez_genes <- as.character(FinalAnnotatedVcf$entrezgene)
entrez_genes <- unique(entrez_genes[!is.na(entrez_genes)])

##  pathway Enrichment
enrichPath = enrichPathway(gene=entrez_genes, readable=T, pvalueCutoff = 0.01)





pdf('Bio_Analysis.pdf')
if( isEmpty(geneSetDf) ){
  grid.text("No annotated gene in the vcf file ",x = 0.5, y=0.5 ,gp=gpar(fontsize=12),just = "center")
  grid.newpage()}

p_mutFreq
p_geneScore
dotplot(enrichPath, showCategory=10)
grid.text("Pathway Enrichment",x = 0.5, y=0.99 ,gp=gpar(fontsize=11),just = "center")
emapplot(enrichPath, showCategory=10)
cnetplot(enrichPath, categorySize="pvalue", showCategory=10)
dev.off()


