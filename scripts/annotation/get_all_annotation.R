library(data.table)
args <- commandArgs(T)
go <- args[1]
kegg <- args[2]
out_file <- args[3]

go <- fread(go)
kegg <- fread(kegg)
dt <- merge(x = go, y = kegg, all = T, by = c("gene_id", "gene_name"))
fwrite(dt, file = out_file, na = "", quote = F, sep = "\t")