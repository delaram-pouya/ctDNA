ource('Codes/1_Functions.R')
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



