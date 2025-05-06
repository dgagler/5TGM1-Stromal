#!/usr/bin/env Rscript

# Code for Endothelial-to-Mesenchymal (EndoMT) analysis of mouse 5TGM1 stromal cells. Methodology derived from Kenswil et al., 2021 (10.1016/j.stem.2021.01.006)

# Load libraries
library(Seurat)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(patchwork)

stromals <- readRDS("./Mouse5TGM1_StromalsOnly_AnnotatedObject.rds")

subtype <- subset(stromals, subset = celltypes_ECsubbed %in% c("SEC", "AEC", "MSCs", "OLCs"))

subtype <- ScaleData(subtype)
subtype <- FindVariableFeatures(subtype, nfeatures = 2000)
subtype <- RunPCA(subtype, npcs = 10)
subtype <- FindNeighbors(subtype, dims = 1:10)
subtype <- FindClusters(subtype, resolution = 0.2)
subtype <- RunUMAP(subtype, dims = 1:10)
DimPlot(subtype)

# Subsetting out cells with Cdh5 expression (as per Kenswil et al., 2021)
cdh5 <- subset(stromals, subset = Cdh5 > 0)
# Reclustering
cdh5 <- ScaleData(cdh5)
cdh5 <- FindVariableFeatures(cdh5, nfeatures = 2000)
cdh5 <- RunPCA(cdh5, npcs = 10)
cdh5 <- FindNeighbors(cdh5, dims = 1:10)
cdh5 <- FindClusters(cdh5, resolution = 0.2)
cdh5 <- RunUMAP(cdh5, dims = 1:10)
DimPlot(cdh5)

# Plot it out
VlnPlot(cdh5, features = c("Cdh5", "Pecam1", "Eng", "Emcn", "Prrx1", "Col1a1", "Snai2", "Twist1", "Twist2", "Zeb2", "Yap1", "Cxcl12"), ncol = 4)
DimPlot(cdh5)

# Based on violin plots above
endoMTcells <- subset(cdh5, subset = seurat_clusters %in% c("4", "6"))
endoMTcellnames <- colnames(endoMTcells)

# Select cells with true endoMT expression
plot <- DimPlot(endoMTcells)
select.cells <- CellSelector(plot = plot)

# Assign to objects
Idents(cdh5) <- "seurat_clusters"
Idents(cdh5, cells = select.cells) <- "EndoMT"
cdh5@meta.data$endoMT <- Idents(cdh5)
DimPlot(cdh5)

endoMT.object <- subset(cdh5, subset = endoMT == "EndoMT")

DimPlot(object = cdh5, cells.highlight = select.cells, cols.highlight = "red", cols = "gray", order = TRUE, label = T) + NoLegend()
ggsave("./figures/Cdh5Pos_EndoMT_Highlighted_UMAP.png", height = 4, width = 5)

saveRDS(cdh5, "./5TGM1_Cdh5Pos_EndoMT_Annotated_FromOG_Object.rds")

# Backtracking Endo-MT cells to original object
Idents(subtype) <- "seurat_clusters"
Idents(subtype, cells = select.cells) <- "Endo-MT"
subtype@meta.data$endoMT <- Idents(subtype)
DimPlot(subtype)
saveRDS(subtype, "./5TGM1_BMEC_MSCOLC_EndoMT_Annotated_FromOG_Object.rds")

cdh5 <- readRDS("./5TGM1_Cdh5Pos_EndoMT_Annotated_FromOG_Object.rds")
