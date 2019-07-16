ource('Codes/1_Functions.R')
Initialize()


FinalAnnotatedVcf <- read.csv('Annaotated_VCF.csv')


## make gene table
geneSet <- FinalAnnotatedVcf$hgnc_symbol
geneSet[is.na(geneSet)|geneSet=='']='Not_Annotated'
geneSetDf <- data.frame(table(geneSet))
rownames(geneSetDf) <- geneSetDf$geneSet
colnames(geneSetDf)[1] <- 'Gene' 


listOfDataSets <- ImportDataSets() 
lapply(listOfDataSets, head)
listOfDataSetGenes <- lapply(listOfDataSets, function(x) toupper(x$Gene))


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
head(merged_table)
write.csv(merged_table, 'GeneInfoTable.csv')




