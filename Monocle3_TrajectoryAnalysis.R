library(Seurat)
library(monocle3)
library(monocle)
library(ggplot2)
library(SeuratWrappers)

options(future.globals.maxSize=1000000000000000) # Set max global size so we don't run out of memory

msclin <- readRDS("/Users/gagled01/morganLab/single-cell/5TGM1_Stromal/objects/5TGM1_MSCOLC_Object.rds")
bmecs <- readRDS("/Users/gagled01/morganLab/single-cell/5TGM1_Stromal/objects/5TGM1_BMEC_Object.rds")

subtype <- msclin
#subtype <- bmecs

cds <- as.cell_data_set(subtype)
cds <- estimate_size_factors(cds)

plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

Idents(subtype) <- "seurat_clusters"
list.cluster <- subtype@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# MSCs
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- subtype@reductions$umap@cell.embeddings

cds <- learn_graph(cds, use_partition = T)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F, 
           group_label_size = 5) + facet_grid(. ~ case.control)

#cds <- order_cells(cds)
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "1"]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F) + facet_grid(. ~ condition) 
#ggsave("/Users/gagled01/morganLab/single-cell/5TG_Mouse_Rerun/BMEC_Monocle_TrajectoryUMAP_Pseudotime_ConditionSplit_Cluster2Root.png", height = 5, width = 8)
