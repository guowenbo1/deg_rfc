#!/usr/bin/env Rscript

library(argparser)
library(tidyverse)

Args <- commandArgs(TRUE)
deg = Args[1]
contrast = Args[2]
outall = Args[3]


deg <- read_tsv(deg)

if (!is.null(all)){
  all  <- deg %>% filter(level %in% c("Increased", "Decreased")) %>% .[[1]]
}
all <- as.data.frame(all)
names(all) <- contrast
write_tsv(all,outall)




# if (!is.null(argv$down)){
#   down <- deg %>% filter(level == "down") %>% .[[1]]
#   write(down, argv$down)
# }


# if (!is.null(argv$up)){
#   up <- deg %>% filter(level == "up") %>% .[[1]]
#   write(up, argv$up)
# }