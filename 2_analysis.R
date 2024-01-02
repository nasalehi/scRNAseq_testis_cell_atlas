library(dplyr)
library(Seurat)
library(ggplot2)

setwd("....")

load("integration.RData")
######################################Dimention Reduction##############################################
DefaultAssay(all.integrated) <- "integrated"
all.integrated <- ScaleData(all.integrated, verbose = FALSE)
all.integrated <- RunPCA(all.integrated, npcs = 50, verbose = FALSE)

pdf("./Results/ElbPCA.pdf")
ElbowPlot(all.integrated,ndims = 50)
dev.off()

all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:35)
dim(all.integrated)
#2000 26642
DataTypes <- Idents(all.integrated)
table(Idents(all.integrated))

current.cluster.ids <-  c("D2", "D7", "fetal", "Spermatogenesis1", "Spermatogenesis2", "Y1", "Y11", "Y13", "Y14", "Y7")
new.cluster.ids <-  c("Infant", "Infant", "Fetal", "Puberty", "Puberty", "Childhood", "Childhood", "Peri-puberty", "Peri-puberty", "Childhood")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)

my_levels <-   c("Fetal", "Infant", "Childhood", "Peri-puberty", "Puberty")
all.integrated@active.ident <- factor(x = all.integrated@active.ident, levels = my_levels)
DataTypes2 <- Idents(all.integrated)

######################################Clustering##############################################
all.integrated <- FindNeighbors(all.integrated, dims = 1:35)
all.integrated <- FindClusters(all.integrated, resolution = 0.3)

current.cluster.ids <-  c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17")
new.cluster.ids <-  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
class<- Idents(all.integrated)
table(Idents(all.integrated))
prop.table(table(Idents(all.integrated)))
clusdata <-table(Idents(all.integrated), DataTypes2)
write.csv(clusdata, file= "./Results/clusdata_integrated2.csv")


pdf("./Results/Umap.pdf", width=10)
Idents(all.integrated) <- DataTypes
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE, 
       cols= c("lightskyblue2","deepskyblue4","darkslategray4","darkseagreen","lightpink1"),
       order = c("Fetal", "Infant", "Childhood", "Peri-puberty","Puberty" ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

Idents(all.integrated) <- class
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 8)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
######################################Markers##############################################
DefaultAssay(all.integrated) <- "RNA"
Idents(all.integrated) <- class
###Sertoli
pdf("./Results/Sertoli.pdf", width=10)
FeaturePlot(all.integrated, features = c("GATA4", "AMH", "SOX9" ),cols=c("lightgrey","#FF9900"))
DotPlot(all.integrated, features = c("GATA4", "AMH", "SOX9" ),col.min = 0,cols=c("lightgrey","#FF9900")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
###Leydig
pdf("./Results/Leydig.pdf", width=10)
FeaturePlot(all.integrated, features = c("ARX","C7","ALDH1A1"),cols=c("lightgrey","#9966cc"))
DotPlot(all.integrated, features = c("ARX","C7","ALDH1A1"),col.min = 0,cols=c("lightgrey","#9966cc")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

###Myoid
pdf("./Results/Myoid.pdf", width=10)
FeaturePlot(all.integrated, features = c("ACTA2","MYH11"),cols=c("lightgrey","#33cccc"))
DotPlot(all.integrated, features = c("ACTA2","MYH11"),col.min = 0,cols=c("lightgrey","#33cccc")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
###Macrophages
pdf("./Results/Macrophages.pdf", width=10)
FeaturePlot(all.integrated, features = c("CD163","CD68","MSR1"),cols=c("lightgrey","#669966"))
DotPlot(all.integrated, features = c("CD163","CD68","MSR1"),col.min = 0,cols=c("lightgrey","#669966")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
###Endothelial
pdf("./Results/Endothelial.pdf", width=10)
FeaturePlot(all.integrated, features = c("VWF","SOX17"),cols=c("lightgrey","#ff0033"))
DotPlot(all.integrated, features = c("VWF","SOX17"),col.min = 0,cols=c("lightgrey","#ff0033")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
###Spermatogonia
pdf("./Results/Spermatogonia.pdf", width=10)
FeaturePlot(all.integrated, features = c("PIWIL4","FGFR3","HMGA1","MAGEA4" ),cols=c("lightgrey","#0066ff"))
DotPlot(all.integrated, features = c("PIWIL4","FGFR3","HMGA1","MAGEA4" ),col.min = 0,cols=c("lightgrey","#0066ff")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
###Spermatocyte
pdf("./Results/Spermatocyte.pdf", width=10)
FeaturePlot(all.integrated, features = c("DMC1", "SYCP3", "PIWIL1","OVOL2"),cols=c("lightgrey","#990033"))
DotPlot(all.integrated, features = c("DMC1", "SYCP3", "PIWIL1","OVOL2"),col.min = 0,cols=c("lightgrey","#990033")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
###Spermatid
pdf("./Results/Spermatid.pdf", width=10)
FeaturePlot(all.integrated, features = c("SYCP3","OVOL2", "TEX29","SUN5","SPEM1","ACR","PGK2"),cols=c("lightgrey","#99cc66"))
DotPlot(all.integrated, features = c("SYCP3","OVOL2", "TEX29","SUN5","SPEM1","ACR","PGK2"),col.min = 0,cols=c("lightgrey","#99cc66")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off() 

Idents(all.integrated) <- class
current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)
new.cluster.ids <- c("leydig", "Sertoli-1", "Endothelial", "Undifferentiated SSC", "Early round SPT", "Round SPT-2",
                     "Elongating SPT", "Differentiating SSC", "Diplotene SPC", "Zygotene SPC", "Sertoli-3", 
                     "Macrophage","Sertoli-2", "Round SPT-1", "Leptotene SPC", "Pachytene SPC", "Myoid-1", "Myoid-2")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
cellTypes<-Idents(all.integrated)
DefaultAssay(all.integrated) <- "RNA"

my_levels <-  c("Sertoli-1","Sertoli-2", "Sertoli-3","leydig","Myoid-1", "Myoid-2", "Macrophage","Endothelial","Undifferentiated SSC", "Differentiating SSC", "Leptotene SPC", "Pachytene SPC", "Zygotene SPC", "Diplotene SPC", "Early round SPT", "Round SPT-1", "Round SPT-2", "Elongating SPT")
# Re-level object@ident
all.integrated@active.ident <- factor(x = all.integrated@active.ident, levels = my_levels)
table(Idents(all.integrated), DataTypes)
cellsdata <-table(Idents(all.integrated), DataTypes)
write.csv(cellsdata, file= "./Results/cellsdata_integrated.csv")

pdf("./Results/UMAP2.pdf", width=10)
Idents(all.integrated) <- cellTypes
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 5)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

######################################percent of datasets cells############################################## 
counts <- readxl::read_excel("Results2/BarPlot.xlsx")
p1 <- as.matrix(counts[-1])
colnames(p1)
rownames(p1) <- counts$...1

data_long <- melt(p1, id.vars = rownames(p1))

pdf("./Results2/BarPlot.pdf", width=10)
ggplot(data_long,                  
       aes(x = Var2,
           y = value,
           fill = Var1)) +
  scale_fill_manual(values= c("Fetal"="lightpink1","Infancy"="darkseagreen","Childhood"="darkslategray4","Peri-puberty"="deepskyblue4", "Adulthood"="lightskyblue2")) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme( axis.text.x=element_text(size=12, face = "italic")) +
  NULL
dev.off()
######################################DEGs##############################################
DefaultAssay(all.integrated) <- "RNA"
Idents(all.integrated) <- class 
all.markers <- FindAllMarkers(all.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./Results/DEG_all_integrated.csv")

data_1 <- readxl::read_excel("./Results/DEGs/summary.xlsx")
pdf("./Results2/DEGs/DEG_numbers.pdf", width=10)
positions <- c("Sertoli-1", "Sertoli-2", "Sertoli-3", "leydig", "Myoid-1", "Myoid-2", "Macrophage", "Endothelial")
  data_1 %>%
  ggplot(aes(x = Cell_Types, y = DEG_No, fill = DEG))+ scale_x_discrete(limits = positions) +
  geom_bar(stat = "identity")+ 
  theme(text = element_text(size=20))+
  coord_flip()
dev.off()

####Enrichment Plot
DEGsData <- data.frame(readxl::read_excel("Results/DEGs/enrichment.xlsx"))
DEGsData$BP <- factor(DEGsData$BP,                                   
                  levels = c("aerobic respiration", "ATP synthesis coupled electron transport","cellular response to steroid hormone stimulus","male gonad development",
                    "very-low-density lipoprotein particle clearance", "positive regulation of cholesterol esterification", "phospholipid efflux",
                    "cellular respiration", "mitochondrial electron transport, cytochrome c to oxygen", "mitochondrial ATP synthesis",
                    "elastic fiber assembly", "extracellular matrix organization", "response to estradiol",
                    "platelet aggregation", "angiogenesis1", "muscle contraction", 
                    "detoxification of copper ion", "negative regulation of growth", "cellular oxidant detoxification", 
                    "immune response","inflammatory response", "antigen processing and presentation",
                   "angiogenesis", "positive regulation of gene expression", "response to hypoxia"))

pdf("./Results2/DEGs/Enrich.pdf",width=9)
ggp <- ggplot(DEGsData, aes(BP, P)) + 
  geom_bar(stat = "identity")+                                      
  coord_flip() + theme(text = element_text(size=15))+
  xlim(rev(as.character(DEGsData$BP)))
ggp  
dev.off()

data_1 <- readxl::read_excel("./Results2/DEGs/DEGSummary.xlsx")
pdf("./Results2/DEGs/DEG_No.pdf", width=10)
positions <- c("Sertoli-1", "Sertoli-2", "Sertoli-3", "leydig", "Myoid-1", "Myoid-2", "Macrophage", "Endothelial")
               
  data_1 %>%
  ggplot(aes(x = Cell_Types, y = DEG_No, fill = DEG))+ scale_x_discrete(limits = positions) +
  geom_bar(stat = "identity")+ 
  theme(text = element_text(size=20))+
  coord_flip()
dev.off()

save.image(file = "analysis.RData")
