library(tidyverse)

method <- snakemake@wildcards[["method"]]
type <- snakemake@params[["type"]]
logfc_thres <- as.numeric(snakemake@params[["logfc"]])
value_thres <- as.numeric(snakemake@params[["value"]])
print(type)

data <- read_tsv(snakemake@input[["de_res"]]) %>%
  rename(gene_id = gene)
gene_id2name <- read_tsv(snakemake@input[["gene_id2name"]])

if (method == "DESeq2") {
  data <- left_join(data, gene_id2name) %>% mutate(pvalue = replace_na(pvalue,1)) %>%
    rename(logFC = log2FoldChange, FDR = padj)
}else {
  data <- left_join(data, gene_id2name) %>% rename(pvalue = PValue)
}


# !!sym(type)
data <- data %>% mutate(level = case_when(logFC > logfc_thres & !!sym(type) < value_thres ~ "Increased",
                                          logFC < -logfc_thres & !!sym(type) < value_thres ~ "Decreased",
                                          abs(logFC) <= logfc_thres | !!sym(type) >= value_thres ~ "nonsignificant"))

head(data$level)
fpkm <- read_tsv(snakemake@input[["fpkm"]])
colnames(fpkm)[1] <- "gene_id"
colnames(fpkm)[2:length(colnames(fpkm))] <- str_c(colnames(fpkm)[2:length(colnames(fpkm))], "(FPKM)")
tpm <- read_tsv(snakemake@input[["tpm"]])
colnames(tpm)[1] <- "gene_id"
colnames(tpm)[2:length(colnames(tpm))] <- str_c(colnames(tpm)[2:length(colnames(tpm))], "(TPM)")

data <- left_join(data, fpkm) %>% left_join(tpm)

write_tsv(data, snakemake@output[[1]])
