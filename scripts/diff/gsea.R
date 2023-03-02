library(clusterProfiler)
library(enrichplot)
library(data.table)
library(ggplot2)

Args <- commandArgs(TRUE)
go = Args[1]
kegg = Args[2]
deg = Args[3]
go_gsea = Args[4]
kegg_gsea = Args[5]
outdir = Args[6]
print(outdir)

dir.create(paste0(outdir,"GO/"))
dir.create(paste0(outdir,"KEGG/"))
go <- fread(go)
gene_express <- read.delim(deg,stringsAsFactors = FALSE)
genelist <- gene_express$logFC

names(genelist) <- gene_express$gene_id
genelist <- genelist[order(genelist,decreasing = TRUE)]

gsea_go <- GSEA(genelist,
	TERM2GENE = go[,c(3,1)],
	TERM2NAME = go[,c(3,6)],
	pvalueCutoff = 1, 
	pAdjustMethod = 'BH',
	minGSSize = 10,
	maxGSSize = 1000)

write.table(gsea_go,go_gsea,sep = '\t',row.names = FALSE,quote = FALSE)

count <- 1
for (i in gsea_go$Description) {
  if(is.na(i)){
    count <- count + 1
    next
  }
  p <- gseaplot2(gsea_go, title = i,geneSetID = count, color = "firebrick", pvalue_table = T)
  i_format <- gsub("/", "|", i)
  ggsave(paste0(outdir,"GO/",count,"_", i_format, ".pdf"), p, height = 9, width = 12)
  print(i_format)
  print(paste0(count, "/", length(gsea_go$Description), " has been saved."))
  count <- count + 1
  if(count > 50) {
    break
  }
}

kegg <- fread(kegg)

gsea_kegg <- GSEA(genelist,
    TERM2GENE = kegg[,c("level3_pathway_id","gene_id")],
    TERM2NAME = kegg[,c("level3_pathway_id","level3_pathway_name")],
    pvalueCutoff = 1, 
    pAdjustMethod = 'BH',
    minGSSize = 10,
    maxGSSize = 1000)

write.table(gsea_kegg,kegg_gsea,sep = '\t',row.names = FALSE,quote = FALSE)

count <- 1
for (i in gsea_kegg$Description) {
  if(is.na(i)){
    count <- count + 1
    next
  }
  p <- gseaplot2(gsea_kegg, title = i,geneSetID = count, color = "firebrick", pvalue_table = T)
  i_format <- gsub("/", "|", i)
  ggsave(paste0(outdir,"KEGG/", count,"_",i_format, ".pdf"), p, height = 9, width = 12)
  print(i_format)
  print(paste0(count, "/", length(gsea_kegg$Description), " has been saved."))
  count <- count + 1
  if(count > 50) {
    break
  }
}
