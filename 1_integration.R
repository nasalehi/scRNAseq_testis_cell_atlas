### integration using Seurat4.3.0 standard integration flow
### outputs .rds file containing the integrated seurat object.

library(dplyr)
library(Seurat)
setwd(".....")

#################################Fetal#####################################
###GSE86146
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
fetal <- CreateSeuratObject(counts = p1.data, project = "fetal", min.cells = 3, min.features = 500)
fetal[["percent.mt"]] <- PercentageFeatureSet(fetal, pattern = "^MT-")
              
pdf("./Results/quality1_vlnplot_F.pdf", width=20)
lapply(c(fetal),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off() 
#################################Infancy#####################################   
###GSE124263           
p3.data <- Read10X(data.dir = "Data/GSE124263/Day2/")
dim(p3.data)
#32738  3634
D2 <- CreateSeuratObject(counts = p3.data, project = "D2", min.cells = 3, min.features = 500)
dim(D2)
#18608  3630
p4.data <- Read10X(data.dir = "Data/GSE124263/Day7/")
dim(p4.data)
#32738  5155
D7 <- CreateSeuratObject(counts = p4.data, project = "D7", min.cells = 3, min.features = 500)
dim(D7)
#18956  5151
remove(p3.data)
remove(p4.data)

D2[["percent.mt"]] <- PercentageFeatureSet(D2, pattern = "^MT-")
D7[["percent.mt"]] <- PercentageFeatureSet(D7, pattern = "^MT-")

pdf("./Results/quality1_vlnplot_D.pdf", width=20)
lapply(c(D2,D7),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

D2 <- subset(D2, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 25)
dim(D2)
#18608  3370
D7 <- subset(D7, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 25)
dim(D7)
#18956  3205

pdf("./Results/quality2_vlnplot_D.pdf", width=20)
lapply(c(D2,D7),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()           

#################################Childhood##################################### 
###GSE120506              
p2.data <-read.delim("Data/GSE120506_infant_combined_UMI.txt")
row.names(p2.data) <- p2.data[,1]
dim(p2.data)
#19123  1341
y1 <- CreateSeuratObject(counts = p2.data, project = "y1", min.cells = 3, min.features = 500)
dim(y1)
#17162  1340
remove(p2.data)
y1[["percent.mt"]] <- PercentageFeatureSet(y1, pattern = "^MT-", assay = 'RNA')

pdf("./Results/quality1_vlnplot_y1.pdf", width=20)
lapply(c(y1),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

y1 <- subset(y1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)
dim(y1)
#17162   1309

pdf("./Results/quality2_vlnplot_y1.pdf", width=20)
lapply(c(y1),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
###GSE134144
p5.data <-read.delim("Data/GSE134144_Pubertal_combined_UMI.txt")
row.names(p5.data) <- p5.data[,1]
dim(p5.data)
#33694 12918
#write.table(colnames(p5.data),"p5.txt")
y7 <- CreateSeuratObject(counts = p5.data[,10951:12918], project = "y7", min.cells = 3, min.features = 500)
dim(y7)
#18898  1960
y11 <- CreateSeuratObject(counts = p5.data[,6775:10950], project = "y11", min.cells = 3, min.features = 500)
dim(y11)
#20517  3920

y7[["percent.mt"]] <- PercentageFeatureSet(y7, pattern = "^MT-", assay = 'RNA')
y11[["percent.mt"]] <- PercentageFeatureSet(y11, pattern = "^MT-", assay = 'RNA')  
              
pdf("./Results/quality1_vlnplot_y2.pdf", width=20)
lapply(c(y7,y11),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()   
              
y7 <- subset(y7, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)
dim(y7)
#18898   1840
y11 <- subset(y11, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)
dim(y11)
#20517   836  
pdf("./Results/quality2_vlnplot_y2.pdf", width=20)
lapply(c(y7,y11),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()              
#################################Peri-Puberty#####################################               
y13 <- CreateSeuratObject(counts = p5.data[,2:4052], project = "y13", min.cells = 3, min.features = 500)
dim(y13)
#22553  3910
y14 <- CreateSeuratObject(counts = p5.data[,4053:6774], project = "y14", min.cells = 3, min.features = 500)
dim(y14)
# 25296  2718
remove(p5.data)

y13[["percent.mt"]] <- PercentageFeatureSet(y13, pattern = "^MT-", assay = 'RNA')
y14[["percent.mt"]] <- PercentageFeatureSet(y14, pattern = "^MT-", assay = 'RNA')

pdf("./Results/quality1_vlnplot_y3.pdf", width=20)
lapply(c(y13,y14),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


y13 <- subset(y13, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
dim(y13)
#22553  2310
y14 <- subset(y14, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
dim(y14)
#25296  2596

pdf("./Results/quality2_vlnplot_y3.pdf", width=20)
lapply(c(y13,y14),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#################################Spermatogenesis#####################################  
###GSE109037
Spermatogenesis1.data <- Read10X(data.dir = "Data/GSE109037/AdultHumanSpermatogenesis-reps1-3")
dim(Spermatogenesis1.data)
#33694  7134
Spermatogenesis1 <- CreateSeuratObject(counts = Spermatogenesis1.data , project = "Spermatogenesis1", min.cells = 3, min.features = 200)
dim(Spermatogenesis1)
#29119  7132
remove(Spermatogenesis1.data)
Spermatogenesis1[["percent.mt"]] <- PercentageFeatureSet(Spermatogenesis1, pattern = "^MT-")
###GSE106487
Spermatogenesis2.data <- read.delim("Data/GSE106487/C1.txt")
dim(Spermatogenesis2.data)
#24153  3059
Spermatogenesis2<- CreateSeuratObject(counts = Spermatogenesis2.data , project = "Spermatogenesis2", min.cells = 3, min.features = 200)
dim(Spermatogenesis2)
#22911  3046
remove(Spermatogenesis2.data)
Spermatogenesis2[["percent.mt"]] <- PercentageFeatureSet(Spermatogenesis2, pattern = "^MT-")

pdf("./Results/quality1_vlnplot_s.pdf", width=20)
lapply(c(Spermatogenesis1,Spermatogenesis2),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

Spermatogenesis1 <- subset(Spermatogenesis1, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 25)
dim(Spermatogenesis1)
#29119  7037
Spermatogenesis2 <- subset(Spermatogenesis2, subset = nFeature_RNA > 800 & nFeature_RNA < 12000)
dim(Spermatogenesis2)
#22911  2971

pdf("./Results/quality2_vlnplot_s.pdf", width=20)
lapply(c(Spermatogenesis1,Spermatogenesis2),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
              
save.image(file = "AllData.RData")
#################################Normalization#####################################
load("AllData.RData")        
all.list <- c(fetal,D2,D7, y1, y7,y11,y13,y14,Spermatogenesis1, Spermatogenesis2)
names(all.list) <- c("fetal","D2","D7", "y1", "y7","y11","y13","y14","Spermatogenesis1", "Spermatogenesis2")
for (i in names(all.list)) {
  all.list[[i]] <- NormalizeData(all.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}

#################################Variable#####################################
for (i in names(all.list)) {
  all.list[[i]] <- FindVariableFeatures(all.list[[i]], selection.method = "vst", nfeatures = 2000)
}
#################################removing#####################################
lapply( c(fetal,D2,D7, y1, y7,y11,y13,y14,Spermatogenesis1, Spermatogenesis2),remove)
#################################INTEGRATION#####################################              
k.filter <- min(sapply(all.list, ncol)) + 1
all.anchors <- FindIntegrationAnchors(object.list=all.list,dims = 1:35, k.filter = k.filter)
remove(all.list)
all.integrated <- IntegrateData(anchorset = all.anchors, dims = 1:35, k.weight=10)
dim(all.integrated )
#2000 26642
remove(all.anchors)
save.image(file = "integration.RData")

