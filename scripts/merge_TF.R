library(tidyverse)
library(data.table)
library(magrittr)

args <- commandArgs(T)
in_file <- args[1]
TF_file <- args[2]
out_file <- args[3]


diff_table <- fread(in_file) %>% select('gene_id')
TF_table <- fread(TF_file)

merged <- merge(diff_table, TF_table)

fwrite(merged, out_file, sep="\t", quote=F)
