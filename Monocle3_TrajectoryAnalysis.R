# Code used for trajectory analysis via Monocle3 package (https://cole-trapnell-lab.github.io/monocle3/)
# Essentially using the methodology outlined in https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

# Load libraries
library(Seurat)
library(monocle3)
library(monocle)
library(ggplot2)
library(SeuratWrappers)
options(future.globals.maxSize=1000000000000000) # Set max global size so we don't run out of memory

# Load objects
msclin <- readRDS("./objects/5TGM1_MSCOLC_Object.rds")
bmecs <- readRDS("./objects/5TGM1_BMEC_Object.rds")

# Set to appropriate group for analysis
subtype <- msclin
#subtype <- bmecs

# Create CDS bobject
cds <- as.cell_data_set(subtype)
cds <- estimate_size_factors(cds)

# Generate partition
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# Set clusters
Idents(subtype) <- "seurat_clusters"
list.cluster <- subtype@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- subtype@reductions$umap@cell.embeddings

# Learn graph
cds <- learn_graph(cds, use_partition = T)

# Order cells using cluster 1 as the root
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "1"]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F) + facet_grid(. ~ condition) 
