source('Codes/1_Functions.R')
Initialize()


FinalAnnotatedVcf <- read.csv('Annaotated_VCF.csv')

entrez_genes <- as.character(FinalAnnotatedVcf$entrezgene)
entrez_genes <- unique(entrez_genes[!is.na(entrez_genes)])


##  pathway Enrichment
pdf('PathwayEnrichment.pdf')
enrichPath = enrichPathway(gene=entrez_genes, readable=T, pvalueCutoff = 0.01)
dotplot(enrichPath, showCategory=10)
grid.text("Pathway Enrichment",x = 0.5, y=0.99 ,gp=gpar(fontsize=11),just = "center")
emapplot(enrichPath, showCategory=10)
cnetplot(enrichPath, categorySize="pvalue", showCategory=10)
dev.off()

