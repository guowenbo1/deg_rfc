library(tidyverse)

eggnog <- read_tsv(snakemake@input[["eggnog"]])

gene2ko <- eggnog %>% select(gene_id, gene_name, KEGG_ko) %>% na.omit() %>%
  separate(KEGG_ko, paste0("X", 1:(max(str_count(.$KEGG_ko,","))+1)), sep = ",") %>%
  gather(key = "X", value = "ko", -gene_id, -gene_name) %>%
  mutate(ko = str_replace(ko,"ko:", "")) %>% select(-X) %>% na.omit() %>% unique()

kegg <- read_tsv(snakemake@input[["kegg"]],na = "-")

kegg_anno <- gene2ko %>% left_join(kegg) %>%
  select(-level1_pathway_id, -level2_pathway_id) %>% 
  mutate(level3_pathway_id = str_c("ko", level3_pathway_id))

write_tsv(kegg_anno, snakemake@output[["kegg_anno"]], na = "-")

ko <- kegg %>% select(ko, ko_name, ko_des, ec, level3_pathway_id,level3_pathway_name) %>%
  mutate(pathway_id = str_c("ko",level3_pathway_id)) %>% 
  mutate(pathway = str_c(pathway_id, ",", level3_pathway_name)) %>% 
  group_by(ko, ko_name, ko_des, ec) %>%
  summarise(pathway_num = n(),pathway = str_c(pathway, collapse = ";"))

kegg_anno_gene2pathway <- left_join(gene2ko, ko)
write_tsv(kegg_anno_gene2pathway, snakemake@output[["gene2pathway"]])

kegg_anno_type2gene <- kegg %>% select(level1_pathway_name, level2_pathway_name, ko) %>%
  left_join(gene2ko) %>% na.omit() %>% 
  select(level1_pathway_name, level2_pathway_name, gene_id, gene_name) %>%
  group_by(level1_pathway_name, level2_pathway_name) %>% 
  summarise(gene_num = n(), gene_id = str_c(gene_id, collapse = ";"),
            gene_name = str_c(gene_name, collapse = ";")) %>% 
  group_by(level1_pathway_name) %>%
  arrange(level1_pathway_name, -gene_num)
kegg_anno_type2gene$level2_pathway_name <- factor(kegg_anno_type2gene$level2_pathway_name,
                                                  levels = kegg_anno_type2gene$level2_pathway_name)
write_tsv(kegg_anno_type2gene, snakemake@output[["type2gene"]])


kegg_anno_pathway2gene <- gene2ko %>% left_join(kegg) %>% na.omit() %>%
  select(level3_pathway_id, level3_pathway_name, level1_pathway_name, 
         level2_pathway_name, gene_id, gene_name) %>%
  mutate(level3_pathway_id = str_c("ko", level3_pathway_id)) %>%
  group_by(level3_pathway_id, level3_pathway_name, level1_pathway_name,level2_pathway_name)%>%
  summarise(gene_num=n(), gene_id = str_c(gene_id, collapse = ";"),
            gene_name = str_c(gene_name, collapse = ";"))
write_tsv(kegg_anno_pathway2gene, snakemake@output[["pathway2gene"]])

ggplot(data = kegg_anno_type2gene) + 
  geom_col(aes(x=level2_pathway_name, y = gene_num, fill = level1_pathway_name)) +
  coord_flip() +   theme_bw() + ylab("") + xlab("") +
  theme(
    panel.grid = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )
ggsave(snakemake@output[["fig"]][1], height=8, width = 12)
ggsave(snakemake@output[["fig"]][2], height=8, width = 12)