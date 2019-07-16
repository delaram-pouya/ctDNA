## This script, gets a vcf file as an input and
## annotates it gencode-v29(hg38) reference genome
## NOTE: MAKE SURE THE COORDINATIONS OF VCF FILE IS BASED ON HG38 REFERENCE GENOME

## input: vcf file
## output: annotated csv file with gene symboles

source('annotate/Codes/1_Functions.R')
Initialize()


## import gencode-v29(hg38) annotation file
AnnotationHg38 <- getCleanedAnnotationFile()

## takes the vcf file as an input 
VCF_INPUT = 'annotate/Data/p2_cfDNA1_snv_lev3.vcf'


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
FinalAnnotatedVcf <- read.csv('Annaotated_VCF.csv')



