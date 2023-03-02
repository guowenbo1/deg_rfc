library(clusterProfiler)
library(tidyverse)
library(viridis)
library(tidyr)
library(topGO)
library(DOSE)

Args <- commandArgs(TRUE)

label = Args[1]
contrast = Args[2]
input_deg = Args[3]
input_go = Args[4]
output_res = Args[5]
out_barplot = Args[6]
out_dotplot = Args[7]
outdir = Args[8]
dir.create(outdir, recursive = TRUE)
# label <- snakemake@wildcards[["label"]]
# contrast <- snakemake@wildcards[["contrast"]]

title <- paste0(contrast, "/", label)
go <- read_tsv(input_go)
go2ont <- go[,c(3,5)] %>% unique()
deg <- read_tsv(input_deg)

if (label == "all") {
  gene_list <- filter(deg, level %in% c("Increased","Decreased")) %>% .[["gene_id"]]
}else if (label == "up") {
  gene_list <- filter(deg, level== "Increased") %>% .[["gene_id"]]
}else if (label == "down") {
  gene_list <- filter(deg, level== "Decreased") %>% .[["gene_id"]]
}

go_enrich <- enricher(gene_list, pvalueCutoff = 1, qvalueCutoff = 1, 
                      TERM2GENE = go[,c(3,1)],TERM2NAME = go[,c(3,6)])


res <- as.data.frame(go_enrich) %>% na.omit() %>% left_join(go2ont, by = c("ID"="GOID"))

write_tsv(res, output_res)

# res_sig <- res %>% group_by(ONTOLOGY) %>% arrange(ONTOLOGY, desc(Count)) %>% 
#   dplyr::slice(1:10)

res_sig <- res %>% group_by(ONTOLOGY) %>% arrange(ONTOLOGY, pvalue) %>% 
  dplyr::slice(1:10)

res_sig$Description <- factor(res_sig$Description, levels = rev(res_sig$Description))
# res_sig$ID <- factor(res_sig$ID, levels = rev(res_sig$ID))

ggplot(res_sig, aes(x= Description, y = Count, fill = -log10(pvalue))) + geom_col(colour="#666666") +
  theme_bw() + coord_flip() + facet_grid(ONTOLOGY~.,scales="free") +   
  xlab("") + ylab("Num of Genes") + ggtitle(title) +
  labs(title= "GO term") +
  theme(panel.grid = element_blank()) + 
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),legend.title=element_text(size=13),axis.text=element_text(size=12,color  = "black")) + 
  scale_fill_gradientn(colours=viridis_pal(option = "D")(20)[5:20]) 

ggsave(out_barplot, width = 16, height = 10)

a<- res_sig %>% separate(GeneRatio, c("n1", "n2"), "/")
a$n1 <- as.numeric(a$n1)
a$n2 <- as.numeric(a$n2)
a$GeneRatio <- a$n1/a$n2
ggplot(a, aes(x= Description, y = GeneRatio)) + 
  geom_point(aes(color = -log10(pvalue), size = Count)) +
  theme_bw() + coord_flip() + 
  facet_grid(ONTOLOGY~.,scales="free") + 
  xlab("") + ylab("GeneRatio") + ggtitle(title) +
  labs(title= "GO term") +
  theme(panel.grid = element_blank()) + geom_point(aes(size = Count),shape=21) +
  theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),legend.title=element_text(size=15),axis.text=element_text(size=14,color  = "black")) + 
  scale_color_gradientn(colours=viridis_pal(option = "D")(20)[5:20]) 

ggsave(out_dotplot, width = 16, height = 10)

#|
# ==============================================================================
# 有向无环图
# ==============================================================================
#|

#重新定义enricher函数(BP)
enricher_internal1<-function (gene, pvalueCutoff, pAdjustMethod = "BH", universe = NULL, 
          minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, USER_DATA) 
{
  
  gene <- as.character(unique(gene))
  qExtID2TermID <- DOSE:::EXTID2TERMID(gene, USER_DATA)
  qTermID <- unlist(qExtID2TermID)
  if (is.null(qTermID)) {
    message("--> No gene can be mapped....")
    p2e <- get("PATHID2EXTID", envir = USER_DATA)
    sg <- unlist(p2e[1:10])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene ID: ", paste0(sg, 
                                                   collapse = ","))
    message("--> return NULL...")
    return(NULL)
  }
  qExtID2TermID.df <- data.frame(extID = rep(names(qExtID2TermID), 
                                             times = lapply(qExtID2TermID, length)), termID = qTermID)
  qExtID2TermID.df <- unique(qExtID2TermID.df)
  qTermID2ExtID <- with(qExtID2TermID.df, split(as.character(extID), 
                                                as.character(termID)))
  extID <- DOSE:::ALLEXTID(USER_DATA)
  if (missing(universe)) 
    universe <- NULL
  if (!is.null(universe)) {
    if (is.character(universe)) {
      extID <- intersect(extID, universe)
    }
    else {
      message("`universe` is not in character and will be ignored...")
    }
  }
  qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)
  qTermID <- unique(names(qTermID2ExtID))
  termID2ExtID <- DOSE:::TERMID2EXTID(qTermID, USER_DATA)
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)
  geneSets <- termID2ExtID
  idx <- DOSE:::get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)
  if (sum(idx) == 0) {
    msg <- paste("No gene set have size >", minGSSize, 
                 "...")
    message(msg)
    message("--> return NULL...")
    return(NULL)
  }
  termID2ExtID <- termID2ExtID[idx]
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))
  k <- sapply(qTermID2ExtID, length)
  k <- k[qTermID]
  M <- sapply(termID2ExtID, length)
  M <- M[qTermID]
  N <- rep(length(extID), length(M))
  n <- rep(length(qExtID2TermID), length(M))
  args.df <- data.frame(numWdrawn = k - 1, numW = M, numB = N - 
                          M, numDrawn = n)
  pvalues <- apply(args.df, 1, function(n) phyper(n[1], n[2], 
                                                  n[3], n[4], lower.tail = FALSE))
  GeneRatio <- apply(data.frame(a = k, b = n), 1, function(x) paste(x[1], 
                                                                    "/", x[2], sep = "", collapse = ""))
  BgRatio <- apply(data.frame(a = M, b = N), 1, function(x) paste(x[1], 
                                                                  "/", x[2], sep = "", collapse = ""))
  Over <- data.frame(ID = as.character(qTermID), GeneRatio = GeneRatio, 
                     BgRatio = BgRatio, pvalue = pvalues, stringsAsFactors = FALSE)
  p.adj <- p.adjust(Over$pvalue, method = pAdjustMethod)
  qobj <- tryCatch(qvalue(p = Over$pvalue, lambda = 0.05, pi0.method = "bootstrap"), 
                   error = function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  }
  else {
    qvalues <- NA
  }
  geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse = "/"))
  geneID <- geneID[qTermID]
  Over <- data.frame(Over, p.adjust = p.adj, qvalue = qvalues, 
                     geneID = geneID, Count = k, stringsAsFactors = FALSE)
  Description <- DOSE:::TERM2NAME(qTermID, USER_DATA)
  if (length(qTermID) != length(Description)) {
    idx <- qTermID %in% names(Description)
    Over <- Over[idx, ]
  }
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1, nc, 2:(nc - 1))]
  Over <- Over[order(pvalues), ]
  Over$ID <- as.character(Over$ID)
  Over$Description <- as.character(Over$Description)
  row.names(Over) <- as.character(Over$ID)
  x <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff, 
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, 
           gene = as.character(gene), universe = extID, geneSets = geneSets, 
           organism = "UNKNOWN", keytype = "UNKNOWN", 
           ontology = "BP", readable = FALSE)
  return(x)
}
enricher1<-function (gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                     universe, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, 
                     TERM2GENE, TERM2NAME = NA) 
{
  USER_DATA <- clusterProfiler:::build_Anno(TERM2GENE, TERM2NAME)
  enricher_internal1(gene = gene, pvalueCutoff = pvalueCutoff, 
                     pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize, 
                     maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = USER_DATA)
}



#重新定义enricher函数(CC)
enricher_internal2<-function (gene, pvalueCutoff, pAdjustMethod = "BH", universe = NULL, 
          minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, USER_DATA) 
{
  
  gene <- as.character(unique(gene))
  qExtID2TermID <- DOSE:::EXTID2TERMID(gene, USER_DATA)
  qTermID <- unlist(qExtID2TermID)
  if (is.null(qTermID)) {
    message("--> No gene can be mapped....")
    p2e <- get("PATHID2EXTID", envir = USER_DATA)
    sg <- unlist(p2e[1:10])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene ID: ", paste0(sg, 
                                                   collapse = ","))
    message("--> return NULL...")
    return(NULL)
  }
  qExtID2TermID.df <- data.frame(extID = rep(names(qExtID2TermID), 
                                             times = lapply(qExtID2TermID, length)), termID = qTermID)
  qExtID2TermID.df <- unique(qExtID2TermID.df)
  qTermID2ExtID <- with(qExtID2TermID.df, split(as.character(extID), 
                                                as.character(termID)))
  extID <- DOSE:::ALLEXTID(USER_DATA)
  if (missing(universe)) 
    universe <- NULL
  if (!is.null(universe)) {
    if (is.character(universe)) {
      extID <- intersect(extID, universe)
    }
    else {
      message("`universe` is not in character and will be ignored...")
    }
  }
  qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)
  qTermID <- unique(names(qTermID2ExtID))
  termID2ExtID <- DOSE:::TERMID2EXTID(qTermID, USER_DATA)
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)
  geneSets <- termID2ExtID
  idx <- DOSE:::get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)
  if (sum(idx) == 0) {
    msg <- paste("No gene set have size >", minGSSize, 
                 "...")
    message(msg)
    message("--> return NULL...")
    return(NULL)
  }
  termID2ExtID <- termID2ExtID[idx]
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))
  k <- sapply(qTermID2ExtID, length)
  k <- k[qTermID]
  M <- sapply(termID2ExtID, length)
  M <- M[qTermID]
  N <- rep(length(extID), length(M))
  n <- rep(length(qExtID2TermID), length(M))
  args.df <- data.frame(numWdrawn = k - 1, numW = M, numB = N - 
                          M, numDrawn = n)
  pvalues <- apply(args.df, 1, function(n) phyper(n[1], n[2], 
                                                  n[3], n[4], lower.tail = FALSE))
  GeneRatio <- apply(data.frame(a = k, b = n), 1, function(x) paste(x[1], 
                                                                    "/", x[2], sep = "", collapse = ""))
  BgRatio <- apply(data.frame(a = M, b = N), 1, function(x) paste(x[1], 
                                                                  "/", x[2], sep = "", collapse = ""))
  Over <- data.frame(ID = as.character(qTermID), GeneRatio = GeneRatio, 
                     BgRatio = BgRatio, pvalue = pvalues, stringsAsFactors = FALSE)
  p.adj <- p.adjust(Over$pvalue, method = pAdjustMethod)
  qobj <- tryCatch(qvalue(p = Over$pvalue, lambda = 0.05, pi0.method = "bootstrap"), 
                   error = function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  }
  else {
    qvalues <- NA
  }
  geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse = "/"))
  geneID <- geneID[qTermID]
  Over <- data.frame(Over, p.adjust = p.adj, qvalue = qvalues, 
                     geneID = geneID, Count = k, stringsAsFactors = FALSE)
  Description <- DOSE:::TERM2NAME(qTermID, USER_DATA)
  if (length(qTermID) != length(Description)) {
    idx <- qTermID %in% names(Description)
    Over <- Over[idx, ]
  }
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1, nc, 2:(nc - 1))]
  Over <- Over[order(pvalues), ]
  Over$ID <- as.character(Over$ID)
  Over$Description <- as.character(Over$Description)
  row.names(Over) <- as.character(Over$ID)
  x <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff, 
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, 
           gene = as.character(gene), universe = extID, geneSets = geneSets, 
           organism = "UNKNOWN", keytype = "UNKNOWN", 
           ontology = "CC", readable = FALSE)
  return(x)
}
enricher2<-function (gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                     universe, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, 
                     TERM2GENE, TERM2NAME = NA) 
{
  USER_DATA <- clusterProfiler:::build_Anno(TERM2GENE, TERM2NAME)
  enricher_internal2(gene = gene, pvalueCutoff = pvalueCutoff, 
                     pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize, 
                     maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = USER_DATA)
}


#重新定义enricher函数(MF)
enricher_internal3<-function (gene, pvalueCutoff, pAdjustMethod = "BH", universe = NULL, 
          minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, USER_DATA) 
{
  
  gene <- as.character(unique(gene))
  qExtID2TermID <- DOSE:::EXTID2TERMID(gene, USER_DATA)
  qTermID <- unlist(qExtID2TermID)
  if (is.null(qTermID)) {
    message("--> No gene can be mapped....")
    p2e <- get("PATHID2EXTID", envir = USER_DATA)
    sg <- unlist(p2e[1:10])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene ID: ", paste0(sg, 
                                                   collapse = ","))
    message("--> return NULL...")
    return(NULL)
  }
  qExtID2TermID.df <- data.frame(extID = rep(names(qExtID2TermID), 
                                             times = lapply(qExtID2TermID, length)), termID = qTermID)
  qExtID2TermID.df <- unique(qExtID2TermID.df)
  qTermID2ExtID <- with(qExtID2TermID.df, split(as.character(extID), 
                                                as.character(termID)))
  extID <- DOSE:::ALLEXTID(USER_DATA)
  if (missing(universe)) 
    universe <- NULL
  if (!is.null(universe)) {
    if (is.character(universe)) {
      extID <- intersect(extID, universe)
    }
    else {
      message("`universe` is not in character and will be ignored...")
    }
  }
  qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)
  qTermID <- unique(names(qTermID2ExtID))
  termID2ExtID <- DOSE:::TERMID2EXTID(qTermID, USER_DATA)
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)
  geneSets <- termID2ExtID
  idx <- DOSE:::get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)
  if (sum(idx) == 0) {
    msg <- paste("No gene set have size >", minGSSize, 
                 "...")
    message(msg)
    message("--> return NULL...")
    return(NULL)
  }
  termID2ExtID <- termID2ExtID[idx]
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))
  k <- sapply(qTermID2ExtID, length)
  k <- k[qTermID]
  M <- sapply(termID2ExtID, length)
  M <- M[qTermID]
  N <- rep(length(extID), length(M))
  n <- rep(length(qExtID2TermID), length(M))
  args.df <- data.frame(numWdrawn = k - 1, numW = M, numB = N - 
                          M, numDrawn = n)
  pvalues <- apply(args.df, 1, function(n) phyper(n[1], n[2], 
                                                  n[3], n[4], lower.tail = FALSE))
  GeneRatio <- apply(data.frame(a = k, b = n), 1, function(x) paste(x[1], 
                                                                    "/", x[2], sep = "", collapse = ""))
  BgRatio <- apply(data.frame(a = M, b = N), 1, function(x) paste(x[1], 
                                                                  "/", x[2], sep = "", collapse = ""))
  Over <- data.frame(ID = as.character(qTermID), GeneRatio = GeneRatio, 
                     BgRatio = BgRatio, pvalue = pvalues, stringsAsFactors = FALSE)
  p.adj <- p.adjust(Over$pvalue, method = pAdjustMethod)
  qobj <- tryCatch(qvalue(p = Over$pvalue, lambda = 0.05, pi0.method = "bootstrap"), 
                   error = function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  }
  else {
    qvalues <- NA
  }
  geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse = "/"))
  geneID <- geneID[qTermID]
  Over <- data.frame(Over, p.adjust = p.adj, qvalue = qvalues, 
                     geneID = geneID, Count = k, stringsAsFactors = FALSE)
  Description <- DOSE:::TERM2NAME(qTermID, USER_DATA)
  if (length(qTermID) != length(Description)) {
    idx <- qTermID %in% names(Description)
    Over <- Over[idx, ]
  }
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1, nc, 2:(nc - 1))]
  Over <- Over[order(pvalues), ]
  Over$ID <- as.character(Over$ID)
  Over$Description <- as.character(Over$Description)
  row.names(Over) <- as.character(Over$ID)
  x <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff, 
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, 
           gene = as.character(gene), universe = extID, geneSets = geneSets, 
           organism = "UNKNOWN", keytype = "UNKNOWN", 
           ontology = "MF", readable = FALSE)
  return(x)
}
enricher3<-function (gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                     universe, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, 
                     TERM2GENE, TERM2NAME = NA) 
{
  USER_DATA <- clusterProfiler:::build_Anno(TERM2GENE, TERM2NAME)
  enricher_internal3(gene = gene, pvalueCutoff = pvalueCutoff, 
                     pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize, 
                     maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = USER_DATA)
}

#使用BP绘制有向无环图
go_BP<- go[which(go$ONTOLOGY=="BP"),]
go_enrich_BP <- enricher1(gene_list, pvalueCutoff = 0.05 , qvalueCutoff = 1, 
                      TERM2GENE = go_BP[,c(3,1)],TERM2NAME = go_BP[,c(3,6)])

file1 = paste(outdir,"DAGplot_BP.pdf", sep="/")
pdf(file1)
plotGOgraph(go_enrich_BP)
dev.off

#使用CC绘制有向无环图
go_CC<- go[which(go$ONTOLOGY=="CC"),]
go_enrich_CC <- enricher2(gene_list, pvalueCutoff = 0.05 , qvalueCutoff = 1, 
                      TERM2GENE = go_CC[,c(3,1)],TERM2NAME = go_CC[,c(3,6)])

file1 = paste(outdir,"DAGplot_CC.pdf", sep="/")
pdf(file1)
plotGOgraph(go_enrich_CC)
dev.off()

#使用MF绘制有向无环图
go_MF<- go[which(go$ONTOLOGY=="MF"),]
go_enrich_MF <- enricher3(gene_list, pvalueCutoff = 0.05 , qvalueCutoff = 1, 
                      TERM2GENE = go_MF[,c(3,1)],TERM2NAME = go_MF[,c(3,6)])

file1 = paste(outdir,"DAGplot_MF.pdf", sep="/")
pdf(file1)
plotGOgraph(go_enrich_MF)
dev.off()
