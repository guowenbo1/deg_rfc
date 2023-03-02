library(dplyr)
library(data.table)
library(tidyverse)
args <- commandArgs(T)

if(length(args)==0)
{
        print("Rscript <gene_count> <gene_fpkm> <gene_tpm> <name_anno> <count_name_anno> <fpkm_name_anno> <tpm_name_anno>")
        quit()
}

gene_count = args[1]
gene_fpkm = args[2]
gene_tpm = args[3]
name_anno = args[4]

count_name_anno = args[5]
fpkm_name_anno = args[6]
tpm_name_anno = args[7]

count <- fread(gene_count)
fpkm <- fread(gene_fpkm)
tpm <- fread(gene_tpm)
nameanno <- fread(name_anno)

names(count)[1] <- names(nameanno)[1]
names(fpkm)[1] <- names(nameanno)[1]
names(tpm)[1] <- names(nameanno)[1]

count_anno <- merge(count,nameanno,all.x=TRUE)
fpkm_anno <- merge(fpkm,nameanno,all.x=TRUE)
tpm_anno <- merge(tpm,nameanno,all.x=TRUE)

count1 <- count_anno %>% select(gene_id,gene_name,everything())
fpkm1 <- fpkm_anno %>% select(gene_id,gene_name,everything())
tpm1 <- tpm_anno %>% select(gene_id,gene_name,everything())

readr::write_tsv(count1 , count_name_anno)
readr::write_tsv(fpkm1 , fpkm_name_anno)
readr::write_tsv(tpm1 , tpm_name_anno)