library(data.table)
library(tidyr)
library(dplyr)
library(readr)

Args <- commandArgs(TRUE)
deg <- Args[1]
diff_anno <- Args[2]

degs <- fread(deg)

gene_list <- degs[level %in% c("Increased","Decreased")]
write_tsv(gene_list,diff_anno)
