library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd('/Users/johnny/Documents/17p_stemness/img_score_validate')


ranks <- read.table("ranked_gene_list_stemness_markers_entrez.rnk",
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$SCORE, ranks$ENTREZID)
str(ranks)

#gene <- names(ranks)
#gene.df <- bitr(gene, fromType = "SYMBOL",
#                toType = c("ENTREZID"),
#                OrgDb = org.Hs.eg.db)

#write.table(gene.df, file = "gene_id_map.txt", row.names = F, sep = "\t", quote = F)



mkk2 <- gseKEGG(geneList = ranks,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
pdf("kegg_enrichment.pdf", width = 10, height = 10)
dotplot(mkk2, showCategory=40)
dev.off()


