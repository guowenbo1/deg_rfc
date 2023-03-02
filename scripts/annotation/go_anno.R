library(GO.db)
library(magrittr)
library(stringr)
library(ggplot2)

getGOLevel <- function(ont, level) {
  switch(ont,
         MF = {
           topNode <- "GO:0003674"
           Children <- GOMFCHILDREN
         },
         BP = {
           topNode <- "GO:0008150"
           Children <- GOBPCHILDREN
         },
         CC = {
           topNode <- "GO:0005575"
           Children <- GOCCCHILDREN
         }
  )
  
  max_level <- max(level)
  if (any(level == 1)) {
    all_nodes <- topNode
  } else {
    all_nodes <- c()
  }
  
  Node <- topNode
  for (i in seq_len(max_level-1)) {
    Node <- mget(Node, Children, ifnotfound=NA)
    Node <- unique(unlist(Node))
    Node <- as.vector(Node)
    Node <- Node[!is.na(Node)]
    if ((i+1) %in% level) {
      all_nodes <- c(all_nodes, Node)
    }
  }
  return(all_nodes)
}

godb <- select(GO.db, keys(GO.db), columns(GO.db))

go_level2 <- union(getGOLevel("MF", 2), getGOLevel("CC", 2)) %>% union(getGOLevel("BP", 2))


eggnog <- readr::read_tsv(snakemake@input[["eggnog"]])

gene2go <- eggnog %>% dplyr::select(gene_id, gene_name, GOs) %>% na.omit() %>%
  tidyr::separate(col = GOs,into=paste0("X", 1:(max(str_count(.$GOs,","))+1)),sep = ",") %>%
  tidyr::gather(key = "X", value = "GOID", -gene_id, -gene_name) %>%
  dplyr::select(gene_id, gene_name, GOID) %>%na.omit() %>% base::unique()

go_annotation <- dplyr::left_join(gene2go, godb)

go_annotation_level2 <- dplyr::filter(gene2go ,GOID %in% go_level2) %>%
  dplyr::left_join(godb)

readr::write_tsv(go_annotation, snakemake@output[["go_anno"]])
readr::write_tsv(go_annotation_level2, snakemake@output[["go_anno_level2"]])


go2gene_level2 <- go_annotation_level2 %>% 
  dplyr::group_by(GOID,DEFINITION,ONTOLOGY,TERM) %>%
  dplyr::summarise(count = length(gene_id),
                   gene_id = paste0(gene_id,collapse = ";"), 
                   gene_name = paste0(gene_name, collapse = ";"))

readr::write_tsv(go2gene_level2, snakemake@output[["go2gene"]])


data <-go2gene_level2 %>% dplyr::ungroup() %>% 
  dplyr::select(GOID, ONTOLOGY, TERM, count) %>%
  dplyr::group_by(ONTOLOGY) %>% dplyr::arrange(ONTOLOGY, -count)

data$TERM <- factor(data$TERM, levels = data$TERM)
genes_num <- dim(eggnog)[1]
# barplot
p <- ggplot(data, aes(x= TERM, y = count, fill = ONTOLOGY)) + geom_col() +
  scale_y_continuous(sec.axis = sec_axis(~.*100/genes_num, name = "Percent of Genes (%)")) +
  theme_bw() + ylab("Number of Genes") + ggtitle("GO level2") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

ggsave(snakemake@output[["fig"]][1], p, height = 8, width = 20)
ggsave(snakemake@output[["fig"]][2], p, height = 8, width = 20)
