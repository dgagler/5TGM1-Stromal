#!/usr/bin/env Rscript

# Script which performs integration on the 5 mouse GEX libraries and then runs a standard Seurat workflow on them

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
options(future.globals.maxSize=1000000000000000)

# Load list of QC filtered libraries
lib.list <- readRDS("./QC_Filtered_Libs.rds")

# Normalize and find variable features for each lib
lib.list <- lapply(X = lib.list, FUN = function(x) {
    x <- NormalizeData(x) 
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
features <- SelectIntegrationFeatures(object.list = lib.list)
lib.list <- lapply(X = lib.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

lib.anchors <- FindIntegrationAnchors(object.list = lib.list, anchor.features = features, reduction = "cca")

# This command creates an 'integrated' data assay
lib.combined <- IntegrateData(anchorset = lib.anchors)

mouseBM <- lib.combined

# Remove individual libs
rm(lib.combined, lib1, lib2, lib3, lib4, lib5, lib6)

DefaultAssay(mouseBM) <- "integrated"

# Standard workflow aka scale data, run PCA (and confirm via elbow plot what ideal is), find neighbors, cluster, and run UMAP
mouseBM <- ScaleData(mouseBM)
mouseBM <- RunPCA(mouseBM, npcs = 30)

ElbowPlot(
  object = mouseBM, 
  ndims = 30
) +
  geom_abline(
    aes(intercept = 1.75, slope = 0, color = "red"),
    show.legend = FALSE
  )

mouseBM <- FindNeighbors(mouseBM, dims = 1:30)
mouseBM <- FindClusters(mouseBM, resolution = 0.2)
mouseBM <- RunUMAP(mouseBM, dims = 1:30)
DimPlot(mouseBM, label = T, group.by = "integrated_snn_res.0.2")

# Save out object
saveRDS(lib.combined, "./5TGM1_IntegratedObject_StandardWorkflow.rds")
