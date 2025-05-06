#!/usr/bin/env Rscript

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load data
lib1.data <- Read10X('./082322-1/filtered_feature_bc_matrix')
lib2.data <- Read10X('./082322-2/filtered_feature_bc_matrix')
lib3.data <- Read10X('./082322-3/filtered_feature_bc_matrix')
lib4.data <- Read10X('./other/filtered_feature_bc_matrix')
lib5.data <- Read10X('./090122-1/filtered_feature_bc_matrix')
lib6.data <- Read10X('./090122-2/filtered_feature_bc_matrix')

# Rename hashtags
lib4.data$`Antibody Capture`@Dimnames[[1]][1:8] <- c("hashtag1", "hashtag2", "hashtag3", "hashtag4", "hashtag5", "hashtag6", "hashtag7", "hashtag8")

# Create Seurat objects
lib1 <- CreateSeuratObject(counts = lib1.data$`Gene Expression`, min.cells=3, min.features=200)
lib2 <- CreateSeuratObject(counts = lib2.data$`Gene Expression`, min.cells=3, min.features=200)
lib3 <- CreateSeuratObject(counts = lib3.data$`Gene Expression`, min.cells=3, min.features=200)
lib4 <- CreateSeuratObject(counts = lib4.data$`Gene Expression`, min.cells=3, min.features=200)
lib5 <- CreateSeuratObject(counts = lib5.data$`Gene Expression`, min.cells=3, min.features=200)
lib6 <- CreateSeuratObject(counts = lib6.data$`Gene Expression`, min.cells=3, min.features=200)

# Fix mismatch in cell count between GEX and ADT
# Lib1
lib1.diff <- setdiff(colnames(lib1.data$`Antibody Capture`), colnames(lib1)) # Getting the unmatched cell
lib1.data$`Antibody Capture` <- lib1.data$`Antibody Capture`[, !colnames(lib1.data$`Antibody Capture`) %in% lib1.diff]
# Lib2
lib2.diff <- setdiff(colnames(lib2.data$`Antibody Capture`), colnames(lib2)) # Getting the unmatched cell
lib2.data$`Antibody Capture` <- lib2.data$`Antibody Capture`[, !colnames(lib2.data$`Antibody Capture`) %in% lib2.diff]
# Lib3
lib3.diff <- setdiff(colnames(lib3.data$`Antibody Capture`), colnames(lib3)) # Getting the unmatched cell
lib3.data$`Antibody Capture` <- lib3.data$`Antibody Capture`[, !colnames(lib3.data$`Antibody Capture`) %in% lib3.diff]
# Lib4
lib4.diff <- setdiff(colnames(lib4.data$`Antibody Capture`), colnames(lib4)) # Getting the unmatched cell
lib4.data$`Antibody Capture` <- lib4.data$`Antibody Capture`[, !colnames(lib4.data$`Antibody Capture`) %in% lib4.diff]
# Lib5
lib5.diff <- setdiff(colnames(lib5.data$`Antibody Capture`), colnames(lib5)) # Getting the unmatched cell
lib5.data$`Antibody Capture` <- lib5.data$`Antibody Capture`[, !colnames(lib5.data$`Antibody Capture`) %in% lib5.diff]
# Lib6
lib6.diff <- setdiff(colnames(lib6.data$`Antibody Capture`), colnames(lib6)) # Getting the unmatched cell
lib6.data$`Antibody Capture` <- lib6.data$`Antibody Capture`[, !colnames(lib6.data$`Antibody Capture`) %in% lib6.diff]

# Create hashtag assay
lib1[["HTO"]] <- CreateAssayObject(counts = lib1.data$`Antibody Capture` + 1, colnames = (x = lib1)) # needed to add 1 to make HTODemux() run
lib2[["HTO"]] <- CreateAssayObject(counts = lib2.data$`Antibody Capture` + 1, colnames = (x = lib2))
lib3[["HTO"]] <- CreateAssayObject(counts = lib3.data$`Antibody Capture` + 1, colnames = (x = lib3))
lib4[["HTO"]] <- CreateAssayObject(counts = lib4.data$`Antibody Capture` + 1, colnames = (x = lib4))
lib5[["HTO"]] <- CreateAssayObject(counts = lib5.data$`Antibody Capture` + 1, colnames = (x = lib5))
lib6[["HTO"]] <- CreateAssayObject(counts = lib6.data$`Antibody Capture` + 1, colnames = (x = lib6))

# Annotate basic metadat
lib1[["library"]] <- "lib1"
lib2[["library"]] <- "lib2"
lib3[["library"]] <- "lib3"
lib4[["library"]] <- "lib4"
lib5[["library"]] <- "lib5"
lib6[["library"]] <- "lib6"
lib1[["case.control"]] <- "control"
lib2[["case.control"]] <- "control"
lib3[["case.control"]] <- "myeloma"
lib4[["case.control"]] <- "myeloma"
lib5[["case.control"]] <- "control"
lib6[["case.control"]] <- "myeloma"

# Add mitochondrial %
# Lib 1
lib1[["percent.mt"]] <- PercentageFeatureSet(lib1, pattern="^mt-")
# Lib 2
lib2[["percent.mt"]] <- PercentageFeatureSet(lib2, pattern="^mt-")
# Lib 3
lib3[["percent.mt"]] <- PercentageFeatureSet(lib3, pattern="^mt-")
# Lib 4
lib4[["percent.mt"]] <- PercentageFeatureSet(lib4, pattern="^mt-")
# Lib 5
lib5[["percent.mt"]] <- PercentageFeatureSet(lib5, pattern="^mt-")
# Lib 6
lib6[["percent.mt"]] <- PercentageFeatureSet(lib6, pattern="^mt-")

# Normalize libraries and demultiplex hashtags
# Lib 1
lib1 <- NormalizeData(object = lib1,
                           assay = "HTO",
                           normalization.method = "CLR")
lib1 <- HTODemux(object = lib1,
                      assay = "HTO",
                      kfunc =  "kmeans",
                      positive.quantile = 0.999)
# Lib 2
lib2 <- NormalizeData(object = lib2,
                           assay = "HTO",
                           normalization.method = "CLR")
lib2 <- HTODemux(object = lib2,
                      assay = "HTO",
                      kfunc =  "kmeans",
                      positive.quantile = 0.999)
# Lib 3
lib3 <- NormalizeData(object = lib3,
                           assay = "HTO",
                           normalization.method = "CLR")
lib3 <- HTODemux(object = lib3,
                      assay = "HTO",
                      kfunc =  "kmeans",
                      positive.quantile = 0.999)
# Lib 4
lib4 <- NormalizeData(object = lib4,
                           assay = "HTO",
                           normalization.method = "CLR")
lib4 <- HTODemux(object = lib4,
                      assay = "HTO",
                      kfunc =  "kmeans",
                      positive.quantile = 0.999)
# Lib 5
lib5 <- NormalizeData(object = lib5,
                           assay = "HTO",
                           normalization.method = "CLR")
lib5 <- HTODemux(object = lib5,
                      assay = "HTO",
                      kfunc =  "kmeans",
                      positive.quantile = 0.999)
# Lib 6
lib6 <- NormalizeData(object = lib6,
                           assay = "HTO",
                           normalization.method = "CLR")
lib6 <- HTODemux(object = lib6,
                      assay = "HTO",
                      kfunc =  "kmeans",
                      positive.quantile = 0.999)

# Get cell counts prior to filtering
prefilter_total <- sum(length(colnames(lib1)), length(colnames(lib2)), length(colnames(lib3)), length(colnames(lib5)), length(colnames(lib6)))

# Filter libraries basd on standard murine threshold of 5% mitochondrial transcripts and also somewhat standard removal of low and high values
# Note that library 4 failed due to low cell count (n = 394) and so is not included
lib1 <- subset(x = lib1,
                        subset = percent.mt < 5.00 & 
                        nCount_RNA > 500 & nFeature_RNA > 500 &
                        nCount_RNA < quantile(nCount_RNA, prob = 1 - 5/100) & 
                        nFeature_RNA < quantile(nFeature_RNA, prob = 1 - 2/100)
                        )
lib2 <- subset(x = lib2,
                        subset = percent.mt < 5.00 & 
                        nCount_RNA > 500 & nFeature_RNA > 500 &
                        nCount_RNA < quantile(nCount_RNA, prob = 1 - 5/100) & 
                        nFeature_RNA < quantile(nFeature_RNA, prob = 1 - 2/100)
                        )
lib3 <- subset(x = lib3,
                        subset = percent.mt < 5.00 & 
                        nCount_RNA > 500 & nFeature_RNA > 500 &
                        nCount_RNA < quantile(nCount_RNA, prob = 1 - 5/100) & 
                        nFeature_RNA < quantile(nFeature_RNA, prob = 1 - 2/100)
                        )
#lib4 <- subset(x = lib4,
#                        subset = percent.mt < 5.00 & 
#                        nCount_RNA > 500 & nFeature_RNA > 500 &
#                        nCount_RNA < quantile(nCount_RNA, prob = 1 - 5/100) & 
#                        nFeature_RNA < quantile(nFeature_RNA, prob = 1 - 2/100)
#                        )

lib5 <- subset(x = lib5,
                        subset = percent.mt < 5.00 & 
                        nCount_RNA > 500 & nFeature_RNA > 500 &
                        nCount_RNA < quantile(nCount_RNA, prob = 1 - 5/100) & 
                        nFeature_RNA < quantile(nFeature_RNA, prob = 1 - 2/100)
                        )

lib6 <- subset(x = lib6,
                        subset = percent.mt < 5.00 & 
                        nCount_RNA > 500 & nFeature_RNA > 500 &
                        nCount_RNA < quantile(nCount_RNA, prob = 1 - 5/100) & 
                        nFeature_RNA < quantile(nFeature_RNA, prob = 1 - 2/100)
                        )

# Get cell count after filtering
postfilter_total <- sum(length(colnames(lib1)), length(colnames(lib2)), length(colnames(lib3)), length(colnames(lib5)), length(colnames(lib6)))
print(paste0("pre-filtering cell counts = ", prefilter_total))
print(paste0("post-filtering cell counts = ", postfilter_total))
print(paste0("fraction lost = ", 1 - postfilter_total/prefilter_total))

# Merge libraries into list
lib.list <- c(lib1, lib2, lib3, lib5, lib6)

# Save out object 
saveRDS(lib.list, "./QC_Filtered_Libs.rds")
