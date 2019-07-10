source('Codes/1_Functions.R')
Initialize()

## takes the input directory and vcf file name as an input 
input_DIR = '~/Desktop/'
vcf_name = 'p2_cfDNA1_snv_lev3.vcf'

## makes a directory to write the outputs in it
outout_DIR = paste0(input_DIR, gsub('.vcf','',vcf_name),'/')
output_prefix = paste0(outout_DIR, gsub('.vcf','',vcf_name))

FinalAnnotatedVcf <- read.csv(paste0(output_prefix,'_annaotatedVCF.csv'))

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


