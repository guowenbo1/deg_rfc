library(data.table)
library(tidyr)
library(dplyr)

args <- commandArgs(T)

deg <- fread(args[1])
diff_table <- deg[level %in% c("Increased","Decreased"),]
tf_table <- fread(args[2])
colnames(tf_table) <- c('gene_id','Family')
out_table <- merge(tf_table, diff_table, on='gene_id', all.y = T) 
head(out_table)
out_table <- out_table %>% select(gene_id, Family) %>% filter(!is.na(Family))

fwrite(out_table, args[3], sep="\t", na="-", quote=F)