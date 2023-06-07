load("analysis.RData")
library(monocle3)
library(htmlwidgets)
nPC <-35
cluster.res <- 0.2
Dim = "2D"
input.dir <- ("./Results/Monocle")
output.dir <- sprintf("./Results/Monocle")
DefaultAssay(all.integrated) <- "integrated"

Idents(all.integrated) <- class
current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)
new.cluster.ids <- c("leydig", "Sertoli-1", "Endothelial", "Undifferentiated SSC", "Early round SPT", "Round SPT-2",
                     "Elongating SPT", "Differentiated SSC", "Diplotene SPC", "Zygotene SPC", "Sertoli-3", 
                     "Macrophage","Sertoli-2", "Round SPT-1", "Leptotene SPC", "Pachytene SPC", "Myoid-1", "Myoid-2")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
cellTypes2<-Idents(all.integrated)

# part one, gene annotations
gene_annotation <- as.data.frame(rownames(all.integrated@reductions[["pca"]]@feature.loadings), row.names = rownames(all.integrated@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
head(gene_annotation)

# part two, cell information
cell_metadata <- as.data.frame(all.integrated@assays[["RNA"]]@counts@Dimnames[[2]], row.names = all.integrated@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
head(cell_metadata)

# part three, counts sparse matrix
New_matrix <- all.integrated@assays[["RNA"]]@data
New_matrix <- New_matrix[rownames(all.integrated@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
New_matrix[1:10,1:10]

### Construct the basic cds object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

### Construct and assign the made up partition
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign the cluster info
list_cluster <- all.integrated@meta.data[[sprintf("ClusterNames_%s_%sPC", cluster.res, nPC)]]
list_cluster <-all.integrated@meta.data$seurat_clusters
list_cluster <-Idents(all.integrated)
names(list_cluster) <- all.integrated@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

### Could be a space-holder, but essentially fills out louvain parameters
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_from_seurat@int_colData@listData$reducedDims$UMAP <-all.integrated@reductions[["umap"]]@cell.embeddings

### Assign feature loading for downstream module analysis
cds_from_seurat@preprocess_aux$gene_loadings <- all.integrated@reductions[["pca"]]@feature.loadings


remove(all.integrated)
### Learn graph, this step usually takes a significant period of time for larger samples
print("Learning graph, which can take a while depends on the sample")
cds_from_seurat <- learn_graph(cds_from_seurat, verbose = FALSE, 
                               learn_graph_control = list(geodesic_distance_ratio=2))


### Plot cluster info with trajectory
print("Plotting clusters")
Dim = "2D"
if (Dim == "2D") {
  pdf(sprintf("%s/clusters.with.trajectory14.%s.pdf", output.dir, Dim), width = 10, height = 10)
  clus <- plot_cells(cds_from_seurat, 
                     color_cells_by = 'cluster',
                     label_groups_by_cluster=TRUE,
                     label_leaves=FALSE,
                     group_label_size=5,
                     label_branch_points=TRUE) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
    theme(legend.position = "right")
  print(clus)
  dev.off()
}

root_cells = choose_cells(cds_from_seurat ,reduction_method ="UMAP", return_list = TRUE)
cds_from_seurat <- order_cells(cds_from_seurat, root_cells = root_cells)

cds_from_seurat_sub <- choose_graph_segments(cds_from_seurat)

pdf("./Results/Monocle/trajectory.pdf", width = 10, height = 10)
plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=1.5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  theme(legend.position = "right")

dev.off()

save.image(file = "monocle.RData")
