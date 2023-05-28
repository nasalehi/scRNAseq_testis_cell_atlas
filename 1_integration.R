### integration using Seurat4.3.0 standard integration flow
### outputs .rds file containing the integrated seurat object.

library(dplyr)
library(Seurat)
setwd(".....")

#################################GSE86146#####################################
filelist <- list.files(path= "Data/GSE86146_RAW/", recursive = TRUE,  pattern = "\\.txt$",full.names = TRUE)
data <- lapply(filelist, read.delim)
dataf[, c(1, 2, 3)]
data[[1]] <- data[[1]][,-2:-51]
all <- Reduce(function(...) merge(..., by="Gene"),data)
dim(all)
#24153  1187
write.table(all, file = "GSE86146_M.txt", sep = "\t",row.names = FALSE)
#######
p1.data <- read.delim("Data/GSE86146_M.txt")
row.names(p1.data) <-p1.data[,1]
dim(p1.data)
#24153  1187
#write.table(colnames(p1.data),"Data/p1.txt")
