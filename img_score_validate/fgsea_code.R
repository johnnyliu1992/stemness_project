library(fgsea)
library(data.table)
library(ggplot2)

data(examplePathways)
data(exampleRanks)
set.seed(42)

setwd('/Users/johnny/Documents/17p_stemness/img_score_validate')

ranks <- read.table("ranked_gene_list_stemness_markers.rnk",
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$rank, ranks$ID)
str(ranks)

pathways <- gmtPathways("stemness_gene_set.gmt")
str(head(pathways))

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

topPathwaysUp <- fgseaRes[padj < 0.1][head(order(NES), n=11), pathway]
#topPathwaysUp <- fgseaRes[head(order(NES), n=11), pathway]


pdf("fgsea_stemness_marker_enrichment.pdf", width = 8, height = 5)

plotGseaTable(pathways[rev(topPathwaysUp)], ranks, fgseaRes, 
              gseaParam=0.5)

dev.off()


