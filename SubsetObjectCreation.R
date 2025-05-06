#!/usr/bin/env Rscript

# Script which generates object subsets of MSC-lineage (MSCs + OLCs) and BMECs (SEC + AEC), reclusters them, and saves them out

# Load libraries
library(Seurat)

# Load stromals
stromals <- readRDS("./Mouse5TGM1_StromalsOnly_AnnotatedObject.rds")
stromals@meta.data$condition <- stromals@meta.data$case.control
stromals@meta.data$condition[stromals@meta.data$condition=="Healthy"] <- "PBS"
stromals@meta.data$condition[stromals@meta.data$condition=="Myeloma"] <- "5TGM1"

# Subset MSCOLCs
mscolcs <- subset(stromals, subset = celltypes_ECsubbed %in% c("MSCs", "OLCs"))
# Recluster
mscolcs <- ScaleData(mscolcs)
mscolcs <- FindVariableFeatures(mscolcs, nfeatures = 2000)
mscolcs <- RunPCA(mscolcs, npcs = 10)
mscolcs <- FindNeighbors(mscolcs, dims = 1:10)
mscolcs <- FindClusters(mscolcs, resolution = 0.2)
mscolcs <- RunUMAP(mscolcs, dims = 1:10)
mscolcs <- subset(mscolcs, subset = seurat_clusters %in% c("0", "1", "2", "3", "4", "5"))
DimPlot(mscolcs)
# Save out
saveRDS(mscolcs, "./objects/5TGM1_MSCOLC_Object.rds")

# Subset BMECs
bmecs <- subset(stromals, subset = celltypes_ECsubbed %in% c("Arterial ECs", "Sinusoidal ECs"))
# Recluster
bmecs <- ScaleData(bmecs)
bmecs <- FindVariableFeatures(bmecs, nfeatures = 2000)
bmecs <- RunPCA(bmecs, npcs = 7)
bmecs <- FindNeighbors(bmecs, dims = 1:7)
bmecs <- FindClusters(bmecs, resolution = 0.15)
bmecs <- RunUMAP(bmecs, dims = 1:7)
DimPlot(bmecs)
# Save out
saveRDS(bmecs, "./objects/5TGM1_BMEC_Object.rds")
