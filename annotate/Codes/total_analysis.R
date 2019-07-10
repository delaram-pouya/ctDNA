## This script, gets a vcf file as an input and
## annotates it gencode-v29(hg38) reference genome
## NOTE: MAKE SURE THE COORDINATIONS OF VCF FILE IS BASED ON HG38 REFERENCE GENOME


source('annotate/Codes/1_Functions.R')
Initialize()

## takes the input directory and vcf file name as an input 
input_DIR = 'annotate/Data/'
vcf_name = 'p2_cfDNA1_snv_lev3.vcf'

## makes a directory to write the outputs in it
outout_DIR = paste0(input_DIR, gsub('.vcf','',vcf_name),'/')
output_prefix = paste0(outout_DIR, gsub('.vcf','',vcf_name))



dir.create(outout_DIR,showWarnings = F)

## import gencode-v29(hg38) annotation file
AnnotationHg38 <- getCleanedAnnotationFile()

## import vcf file
InputVcf <- read.vcfR( paste0(input_DIR, vcf_name), verbose = FALSE )
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



write.csv(FinalAnnotatedVcf, file = paste0(output_prefix,'_annaotatedVCF.csv'), row.names = F)
#FinalAnnotatedVcf <- read.csv(paste0(output_prefix,'_annaotatedVCF.csv'))


######## Score genes

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




##### Gene annotation table

## make gene table
geneSet <- FinalAnnotatedVcf$hgnc_symbol
geneSet[is.na(geneSet)|geneSet=='']='NotAnnot'
geneSetDf <- data.frame(table(geneSet))
rownames(geneSetDf) <- geneSetDf$geneSet
colnames(geneSetDf)[1] <- 'Gene' 


listOfDataSets <- ImportDataSets() 
lapply(listOfDataSets, head)
listOfDataSetGenes <- lapply(listOfDataSets, function(x) toupper(x$Gene))


## preprocess CGC
CGC <- listOfDataSets[['CGC']]
CGC <- subset(CGC, Gene %in% geneSet)
CGC <- subset(CGC, select=c(Gene, Name, Tier, Chr.Band, Tumour.Types.Germline. ,Tumour.Types.Somatic., 
                            Germline,Somatic, Role.in.Cancer, Mutation.Types ,Translocation.Partner))
colnames(CGC)[2:ncol(CGC)] <- paste0(colnames(CGC)[2:ncol(CGC)] , '(CGC)')
head(CGC)


## PanCanSig
PanCanSig <- listOfDataSets$PanCanSig
PanCanSig <- subset(PanCanSig, Gene %in% geneSet)
colnames(PanCanSig)[2] <- 'signaling_pathway(PanCan)' 
head(PanCanSig)    


## PanCan drivers
PanCanDriver <-listOfDataSets$PanCanDriver
colnames(PanCanDriver)[3] <- 'tsg_onco_20.20pred'
PanCanDriver <- subset(PanCanDriver, Gene %in% geneSet)
PanCanDriver <- subset(PanCanDriver, select=c(Gene, Cancer, tsg_onco_20.20pred,Pancan.Frequency))
colnames(PanCanDriver)[2:ncol(PanCanDriver)] <- paste0(colnames(PanCanDriver)[2:ncol(PanCanDriver)] , '(PanCanDriver)')
head(PanCanDriver)

merged_table = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Gene", all.x = TRUE),list(geneSetDf,PanCanDriver,PanCanSig,CGC))
merged_table <- merged_table[order(merged_table$Freq, decreasing = T),]

head(merged_table)
write.csv(merged_table, paste0(output_prefix, '_GeneInfoTable.csv'))



## DGIdb 
DGIdb <- listOfDataSets$DGIdb
DGIdb <- subset(DGIdb, Gene %in% geneSet)
DGIdb <- subset(DGIdb, select=c(Gene, drug_name, drug_chembl_id))
colnames(DGIdb)[2:ncol(DGIdb)] <- paste0(colnames(DGIdb)[2:ncol(DGIdb)] , '(DGIdb)')
head(DGIdb)


# ToDo:
## how to add drug bank drug-ids
## any other columns to add? >> oncoKB ??
## check the dataset websites or papers
## make it as a pipline



####### Enrichments

entrez_genes <- as.character(FinalAnnotatedVcf$entrezgene)
entrez_genes <- unique(entrez_genes[!is.na(entrez_genes)])


## gene Ontology profile
pdf(paste0(output_prefix,'_GeneOntology_profile.pdf'))
ggo = as.data.frame(groupGO(gene = entrez_genes, OrgDb= org.Hs.eg.db, ont= "CC", level= 3, readable = TRUE))
ggo_sub <- subset(ggo, select=c(Description,Count))
ggo_sub <- subset(ggo_sub,Count>0)
ggo_sub$ratio <- round(ggo_sub$Count/length(entrez_genes),3)
ggplot(subset(ggo_sub,ratio>mean(ratio)), aes(x=Description, y=ratio))+ggtitle('gene-ontology')+
  geom_bar(stat="identity",color="black", fill='pink2',alpha=0.8)+coord_flip()+theme_bw()+xlab('')
dev.off()



##  pathway Enrichment
pdf(paste0(output_prefix,'_PathwayEnrichment.pdf'))
enrichPath = enrichPathway(gene=entrez_genes, readable=T, pvalueCutoff = 0.01)
barplot(enrichPath, showCategory=10)
grid.text("Pathway Enrichment",x = 0.5, y=0.99 ,gp=gpar(fontsize=11),just = "center")
dotplot(enrichPath, showCategory=10)
emapplot(enrichPath, showCategory=10)
cnetplot(enrichPath, categorySize="pvalue", showCategory=10)
dev.off()



### Kegg Analysis
pdf(paste0(output_prefix,'_KeggAnalysis.pdf'))
kk = enrichKEGG(gene= entrez_genes, organism= 'hsa',pvalueCutoff = 0.01 ,pAdjustMethod="BH" )
barplot(kk, showCategory=10)
grid.text("Kegg Analysis",x = 0.5, y=0.99 ,gp=gpar(fontsize=11),just = "center")
dotplot(kk, showCategory=10)
emapplot(kk, showCategory=10)
cnetplot(kk, categorySize="genenum",showCategory=10 )
dev.off()



##### GeneOntology 
pdf(paste0(output_prefix,'_GO_EnrichmentAnalysis.pdf'))
ego = enrichGO(gene=entrez_genes,OrgDb= org.Hs.eg.db,ont= "CC",pAdjustMethod = "BH",
               pvalueCutoff= 0.01,readable= TRUE)
barplot(ego, showCategory=10)
grid.text("Gene-Ontology Enrichment Analysis",x = 0.5, y=0.99 ,gp=gpar(fontsize=11),just = "center")
dotplot(ego, showCategory=10)
emapplot(ego, showCategory=10)
cnetplot(ego, categorySize="pvalue", foldChange=entrez_genes,showCategory=10)
dev.off()



