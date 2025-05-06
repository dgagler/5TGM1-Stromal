#!/usr/bin/env Rscript

# Code for manually annotated stromal populations in mouse bone marrow microenvironment. Annotation schema based on Barayawno et al., 2019

# Load libraries
library(Seurat)
library(ggplot2)

# Load all cells
mouseBM <- readRDS("./5TGM1_AllCells_Integrated.rds")
DefaultAssay(mouseBM) <- "RNA"

DimPlot(mouseBM, label = T, group.by = "seurat_clusters")

# Recluster
DefaultAssay(mouseBM) <- "integrated"
mouseBM <- FindNeighbors(mouseBM, dims = 1:30)
mouseBM <- FindClusters(mouseBM, resolution = 0.4)
mouseBM <- RunUMAP(mouseBM, dims = 1:30)
DimPlot(mouseBM, group.by = "seurat_clusters", label= T)

# Mesenchymal stem cells (MSCs)
FeaturePlot(mouseBM, features = "Lepr")
FeaturePlot(mouseBM, features = "Adipoq")
FeaturePlot(mouseBM, features = "Cxcl12")

# Osteo-lineage cells (OLCs)
FeaturePlot(mouseBM, features = "Bglap")
FeaturePlot(mouseBM, features = "Spp1")
FeaturePlot(mouseBM, features = "Sp7")

# Chondrocytes 
FeaturePlot(mouseBM, features = "Col2a1")
FeaturePlot(mouseBM, features = "Sox9")
FeaturePlot(mouseBM, features = "Acan")

# Pericytes
FeaturePlot(mouseBM, features = "Acta2")
FeaturePlot(mouseBM, features = "Myh11")
FeaturePlot(mouseBM, features = "Mcam")

# Endothelial Cells (BMECs)
FeaturePlot(mouseBM, features = "Cdh5")
FeaturePlot(mouseBM, features = "Kdr")
FeaturePlot(mouseBM, features = "Emcn")

# Fibroblasts
FeaturePlot(mouseBM, features = "S100a4")
FeaturePlot(mouseBM, features = "Fn1")
FeaturePlot(mouseBM, features = "Dcn")

# B cells
bcells <- FeaturePlot(mouseBM, features = "Cd79a")

# HSC/GMP
hsc.gmp <- FeaturePlot(mouseBM, features = c("Elane", "Mpo", "S100a9", "Ngp"), ncol = 2)

# Monocyte/macrophages
momacs <- FeaturePlot(mouseBM, features = c("Ccl3", "Tyrobp"))

# Platelets
platelets <- FeaturePlot(mouseBM, features = "Ppbp")

# RBCs
rbcs <- FeaturePlot(mouseBM, features = "Hba-a1")

# Cycling/proliferating cells
cycling <- FeaturePlot(mouseBM, features = c("Car1", "Cdk6", "Peak1", "Mki67"), ncol = 2)

# Oligodendrocytes
oligodendro <- FeaturePlot(mouseBM, features = c("Cdh19", "Plp1"))

# Motor neurons
motor.neurons <- FeaturePlot(mouseBM, features = c("Chodl", "Hs6st3"))

# T cells (none found)
FeaturePlot(mouseBM, features = "Cd4")
FeaturePlot(mouseBM, features = "Nkg7")

# Plasma cells (none found)
FeaturePlot(mouseBM, features = "Ly6k")
FeaturePlot(mouseBM, features = "Slamf7")

# Update annotations
Idents(mouseBM) <- mouseBM@meta.data$seurat_clusters

# NOTE that this requires manual inspection. Adjust these accordingly based on clustering and expression of above markers
mouseBM <- RenameIdents(mouseBM, "0" = "B cells")
mouseBM <- RenameIdents(mouseBM, "1" = "MoMacs")
mouseBM <- RenameIdents(mouseBM, "2" = "B cells")
mouseBM <- RenameIdents(mouseBM, "3" = "Cycling/proliferating")
mouseBM <- RenameIdents(mouseBM, "4" = "Fibroblasts")
mouseBM <- RenameIdents(mouseBM, "5" = "MSC")
mouseBM <- RenameIdents(mouseBM, "6" = "Platelets")
mouseBM <- RenameIdents(mouseBM, "7" = "SEC")
mouseBM <- RenameIdents(mouseBM, "8" = "AEC")
mouseBM <- RenameIdents(mouseBM, "9" = "B cells")
mouseBM <- RenameIdents(mouseBM, "10" = "Pericytes")
mouseBM <- RenameIdents(mouseBM, "11" = "RBC")
mouseBM <- RenameIdents(mouseBM, "12" = "AEC")
mouseBM <- RenameIdents(mouseBM, "13" = "OLC")
mouseBM <- RenameIdents(mouseBM, "14" = "B cells")
mouseBM <- RenameIdents(mouseBM, "15" = "HSC/GMP")
mouseBM <- RenameIdents(mouseBM, "16" = "HSC/GMP")
mouseBM <- RenameIdents(mouseBM, "17" = "RBC")
mouseBM <- RenameIdents(mouseBM, "18" = "Neuronal")
mouseBM <- RenameIdents(mouseBM, "19" = "Oligodendrocytes")

mouseBM@meta.data$celltypes <- Idents(mouseBM)
DimPlot(mouseBM, group.by = "celltypes", label = T)

# Subset stromals
stromals <- subset(mouseBM, subset = celltypes %in% c("Chondrocytes", "OLC", "MSC", "AEC", "SEC", "Fibroblasts", "Pericytes"))
stromals@meta.data$celltypes <- droplevels(stromals@meta.data$celltypes)

# Recluster stromals
stromals <- ScaleData(stromals)
stromals <- FindVariableFeatures(stromals, nfeatures = 2000)
stromals <- RunPCA(stromals, npcs = 25)
ElbowPlot(
  object = stromals, 
  ndims = 25
) +
  geom_abline(
    aes(intercept = 1.75, slope = 0, color = "red"),
    show.legend = FALSE
  )

# Determine percent of variation associated with each PC
pct <- stromals[["pca"]]@stdev / sum(stromals[["pca"]]@stdev) * 100
cumu <- cumsum(pct) # Calculate cumulative percents for each PC

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

if (co1 < co2) {
  n_pc <- co1
  print(paste0("the elbow PC is ", co1))
} else {
  n_pc <- 2
  print(paste0("the elbow point is ", co2))
}

stromals <- FindNeighbors(stromals, dims = 1:21)
stromals <- FindClusters(stromals, resolution = 0.12)
stromals <- RunUMAP(stromals, dims = 1:21)

# Save out
saveRDS(stromals, "./5TGM1_StromalsOnly_Annotated_Object.rds")
