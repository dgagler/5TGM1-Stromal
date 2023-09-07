# Single-cell RNA Analysis of the Bone Marrow Stromal Compartment

**Multiple Myeloma (MM)** is a plasma cell malignancy characterized by abnormal antibody production. After receiving driver mutations in the germinal center, the founding myeloma cell travels to the bone marrow where it alters its local cellular environment to be more conducive to its long-term residency and proliferation. This is termed the myeloma niche. Healthy stem cell, which give rise to immune cells in normal contexts, also exist in privileged bone marrow niches. **Stromal cells are an essential part of both healthy stem cell and myeloma cell niches**, providing structural support, growth and survival signals, and regulating immune cell accesse. While it is appreciated that stromal cells support myeloma progression, the specific ways in which the stromal compartment responds to the myeloma cell is not well understood. **This project aims to better understand the stromal compartment in the context of multiple myeloma.**

# Quick Experimental Overview:
<img width="826" alt="GitPic1" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/f77999e2-8b97-41ee-af38-51d51e703d16">
To study the myeloma niche in a controlled fashion we utilized a mouse model. The 5TGM1 mouse MM model reliably recapitulates both the genetic and physiological features typical of human MM. For this study, 4 mice were injected with 5TGM1 and 4 mice were injected with phosphate buffer saline (PBS) as healthy controls. 1 week after the detection of M-spike proteins in the blood of the 5TGM1 mice the mice were sacrificed and stromal cells were isolated from their bone marrows via flow cytometry and fluoresence-activated cell sorting (FACS). Finally, these cells underwent single-cell RNA sequencing via 10X Genomics. Data was analyzed in R, primarily following a Seurat framework.

# Quality Control
The data gets to me in the form of Cellranger outputs. Cellranger is a software tool which aligns single-cell RNA sequence reads to a reference genome to generate feature-barcode matrices (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). In the context of single-cell RNA data, features represent genes and barcodes represent individual cells. Notably, this experiment also involved sequencing hashtag oligos (HTO), which can be used to identify which sample a particular cell originated from.

For this experiment, our cells were divided across 6 libraries, or sequencing runs. 1 of the 6 libraries had very low cell counts and was immediately excluded from downstream analysis. For each of the 5 remaining libraries, we read in the Cellranger outputs, create a Seurat object from the RNA assay, add in the HTO assay, and add some basic metadata, as shown below:
```{r}
library(Seurat)
library(dplyr)
library(ggplot2)

# Creating a Seurat object based on gene expression data
lib1.data  <- Read10X('/Users/dgagler/5TGM1/lib1/filtered_feature_bc_matrix')
lib1.data <- CreateSeuratObject(counts = lib1.data$`Gene Expression`, min.cells=3, min.features=200) 

# Due to filtering some cells with the min.cells and min.features parameters above, we will have mismatched matrix sizes between our RNA and HTO assays, which we correct for by removing the filtered cells from the HTO assay.

# Getting the unmatched cells
lib1.diff <- setdiff(colnames(lib1.data$`Antibody Capture`), colnames(lib1))
# Removing them
lib1.data$`Antibody Capture` <- lib1.data$`Antibody Capture`[, !colnames(lib1.data$`Antibody Capture`) %in% lib1.diff]
# Adding +1 to the matrix is a quick fix for a downstream error when demultiplexing the cells.
lib1[["HTO"]] <- CreateAssayObject(counts = lib1.data$`Antibody Capture` + 1, colnames = (x = lib1))

# Adding basic metadata
lib1[["library"]] <- "lib1"
lib1[["condition"]] <- "control"

# Computing the percentage of mitochondrial transcripts per cell
# A high percentage of mitochondrial transcripts is often reflective of dead/dying/stressed cells
lib1[["percent.mt"]] <- PercentageFeatureSet(lib1, pattern="^mt-")
```
After preprocessing all 5 libraries, it's time to perform quality filtering. First, we visualize some basic QC metrics including the number of unique UMIs (transcripts) per cell, the number of unique genes per cell, and the % of mitochondrial transcripts. We assess these UMI and gene metrics because having too few UMIs or genes is suggestive of dead or low quality cells and having too many is suggestive of doublets or multiplets (multiple cells in the same sequencing droplet). Quality filtering practices vary, but for this project we elected to exclude cells with:
* greater than 5% mitochondrial transcripts
* less than 500 unique genes or transcripts
* within the top 2% of unique genes
* within the top 5% of unique transcripts
  
<img width="979" alt="GitPic2" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/5c02c9bb-ac82-44e6-b7ad-40046a5783cb">

This filtering left us with 45,030, reflecting a 14.9% loss.

## HTO Data
After QC filtering, we look towards the hashtags. Hashtag oligos (HTOs) are nucleotide barcodes that are uniquely associated with a specific sample. For each library, we normalize the HTO data and then demultiplex the hashtags (determining which cell came from which sample) via Seurat's HTODemux() function.
```{r}
lib1 <- NormalizeData(object = lib1,
                           assay = "HTO",
                           normalization.method = "CLR")
lib1 <- HTODemux(object = lib1,
                      assay = "HTO",
                      kfunc =  "kmeans",
                      positive.quantile = 0.999)
```
Looking at hashtag expression and the number of singlets identified by the demultiplexing, it is obvious that the hashtag data did not perform very well. In particular, these heatmaps show low HTO expression across libraries and the pie charts show the breakdown of singlets (what we want), doublets, and negatives. Given the low quality of the HTO data, we will exclude this information from downstream analysis.

<img width="519" alt="gitpic3" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/40ac29da-188d-4e43-b141-7baf117cd3f1">
<img width="517" alt="gitpic4" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/d85c1291-5286-4761-8842-2736ea1d07ca">

# Integration
Now that our 5 libraries have been quality filtered, it's time to integrate them together. Due to the inherent nature of single-cell library preparation and sequencing technologies, some amount of variation in gene expression is expected between libraries which has nothing to do with the underlying biology. This is commonly termed *batch effect*. There are many ways to correct for this, but here we utilized a Seurat integration method which simultaneously merges together separate data objects and corrects batch effect.

In short, this involves normalizing and finding variable genes for each library, selecting integration features that are highly variable across libraries, scaling and running a PCA on each library using these shared variable features, finding integration anchors within the shared features, and finally integrating the data together using the integration anchors. This will generate a new "integrated" assay in the Seurat object.
```{r}
# list of QC filtered libraries
lib.list <- c(lib1, lib2, lib3, lib4, lib5)

# normalize and find variable features for each library
lib.list <- lapply(X = lib.list, FUN = function(x) {
    x <- NormalizeData(x) 
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = lib.list)
lib.list <- lapply(X = lib.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# find anchors
lib.anchors <- FindIntegrationAnchors(object.list = lib.list, anchor.features = features, reduction = "cca")

# perform integration
integrated.bm <- IntegrateData(anchorset = lib.anchors)
```
# Standard Workflow
After integration, we perform a "standard workflow" on the data. This involves scaling the data by a factor of 10,000, performing a principal component analysis, clustering, and generating a UMAP to visualize our dataset.
```{r}
DefaultAssay(integrated.bm) <- "integrated"

mouseBM <- ScaleData(integrated.bm)
mouseBM <- RunPCA(integrated.bm, npcs = 30)
mouseBM <- FindNeighbors(integrated.bm, dims = 1:30)
mouseBM <- FindClusters(integrated.bm, resolution = 0.2)
mouseBM <- RunUMAP(integrated.bm, dims = 1:30)
DimPlot(integrated.bm, label = T, group.by = "integrated_snn_res.0.2")
```
# Cellular Annotation
We annotated our clusters manually, using apriori knowledge about immune and stromal marker genes. About 2/3 of our 45,030 cells ended up being immune cells, which isn't entirely unexpected as the stromal compartment is a low abundance population, meaning that we had to start with a large number of cells, and our flow sorting methods are not perfect. We removed immune cells from the downstream analysis, resulting in 14,219 high quality stromal cells.
<img width="886" alt="gitpic5" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/72d12785-531c-4f2e-9107-554524cfb11f">

To annotate the stromal populations, we used a set of marker genes derived from Barawyano et al., 2019. Using a combination of violin plots and feature plots, we can assess which clusters are expressing which marker genes. We used *Lepr* expression to identify mesenchymal stem cells (MSC), *Bglap* for osteo-lineage cells (OLC), *Fn1* for fibroblasts, *Sox9* for chondrocytes, *Acta2* for pericytes, and *Cdh5* for endothelial cells.
<img width="956" alt="gitpic6" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/e2862ee7-6290-4b34-ae35-c90e921743f9">

# Bootstrapping Relative Abundance
We are interested in understanding changes in the bone marrow stromal cell populations in response to MM. One way in which these populations differ is in quantity. As an aside, assessing abundances is epistemologically tricky due to the compositional nature of NGS data, but this is not the place to discuss that. See this discussion by Thomas Quinn for more information on the subject - https://academic.oup.com/bioinformatics/article/34/16/2870/4956011

For the time being, however, we will be satisfied with using relative abundances. From the below barplots, you can see that the myeloma condition has relatively more MSC, SEC, and chondrocytes. The chondrocytes are a very small population (0.4-1.2%), however, so we are more interested in the MSC and SEC cells.

<img width="572" alt="gitpic8" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/89074ede-8ec4-4669-8fab-37255376f042">

But we will want to support any findings we may have statistically. This is made more complicated by the fact that our hashtagging data failed. If the hashtagging data had succeeded, we would have 4 replicates in each of our experimental conditions. Since it didn't, we don't. As a workaround, we will employ a common statistical approach: **the bootstrap**.

Bootstrapping is a resampling method that allows for the calculation of standard errors, confidence intervals, p-values, and hypothesis testing. 

In our case, we are interested in knowing whether the observed relative abundance of our stromal cell types are significantly different between 5TGM1 and control mice. To do this, we subset our dataset into 5TGM1 and control, randomly sample (with replacement) from each of our dataset 100x, compute mean relative abundances in our bootstrapped data, and then perform a two-sample T-test to assess the mean relative abundances per cell type.

```{r}
# subset by condition
myeloma <- subset(stromals, subset = condition %in% "myeloma")
healthy <- subset(stromals, subset = condition %in% "healthy")

# get abundances of cell types for resampling
myeloma.abundances <- myeloma@meta.data$celltypes
healthy.abundances <- healthy@meta.data$celltypes

# total cells by condition
n_cells_myeloma <- ncol(myeloma)
n_cells_healthy <- ncol(healthy)

# randomly subsample with replacement and divide by total cells to get relative abundances
myeloma.resamples <- lapply(1:100, function(i) table(sample(myeloma.abundances, replace = T, prob = ))/n_cells_myeloma)
healthy.resamples <- lapply(1:100, function(i) table(sample(healthy.abundances, replace = T, prob = ))/n_cells_healthy)

# parse out relative abundance distributions for each celltype -- MM
mm.resamples.sec <- sapply(myeloma.resamples,"[[",1)
mm.resamples.msc <- sapply(myeloma.resamples,"[[",2)

# parse out relative abundance distributions for each celltype -- Healthy
healthy.resamples.sec <- sapply(healthy.resamples,"[[",1)
healthy.resamples.msc <- sapply(healthy.resamples,"[[",2)
```
Once we've performed our bootstrapping, we'll compute p-values and visualize our data via boxplots using rstatix and patchwork

***note: for the sake of brevity, I've only shown code in the below block for 2 of the 7 cell types
```{r}
library(rstatix)
library(patchwork)

# create dfs for plotting
sec.plotting.df <- melt(data.frame(Myeloma = mm.resamples.sec, Healthy = healthy.resamples.sec, Celltype = "SEC"))
colnames(sec.plotting.df)[2] <- "Condition"

msc.plotting.df <- melt(data.frame(Myeloma = mm.resamples.msc, Healthy = healthy.resamples.msc, Celltype = "MSC"))
colnames(msc.plotting.df)[2] <- "Condition"

# generate individual box plots by cell type
sec.bxp <- ggboxplot(sec.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE) +
  ylab("") + xlab("") + ggtitle("SEC") +
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + NoLegend()
sec.stat.test <- sec.stat.test %>% add_xy_position("Condition")
sec.bxp <- sec.bxp + 
  stat_pvalue_manual(sec.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

msc.bxp <- ggboxplot(msc.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE) +
  ylab("Relative Abundance") + xlab("") + ggtitle("MSC") +
  theme_bw(base_size = 16) +
theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + NoLegend()
msc.stat.test <- msc.stat.test %>% add_xy_position("Condition")
msc.bxp <- msc.bxp + 
  stat_pvalue_manual(msc.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

msc.bxp | olc.bxp | sec.bxp | aec.bxp | fibro.bxp | chondro.bxp | peri.bxp
```
![gitpic7](https://github.com/dgagler/5TGM1-Stromal/assets/31450828/e7bfbea9-c816-4a81-935f-a49cf174d353)

From the above, we can see that the bootstrapping has allowed us to derive p-values for our mean comparisons, which all turned out to be significant. Now we will move forward into more analyses on the populations that are more abundant in the myeloma condition.

# Subcluster Analysis
First, we isolate out our populations of interest. Those we will define as **MSC-lineage** cells, comprised of both our MSC and OLC populations, and **bone marrow endothelial cells (BMEC)**. The reason we include OLC in our MSC-lineage cells is because it is understood that mesenchymal stem cells have the capacity to differentiate into adipocytes (fat cells) and osteo-lineage cells (bone cells), so it makes sense to include all the cells within that known differentiation trajectory.

So, we will subset out our MSC-lineage and BMEC subsets and rerun standard workflows on them. This will allow us to see these populations in higher resolution and shake out subpopulations that may have been concealed in our original clustering analysis which included all stromal cells.

```{r}
msclin <- subset(stromals, subset = celltypes %in% c("MSC", "OLC"))
bmec <- subset(stromals, subset = celltypes %in% c("AEC", "SEC"))

# note - number of PCs and clustering resolution were determined after some trial and error
mscolcs <- ScaleData(msclin)
mscolcs <- FindVariableFeatures(msclin, nfeatures = 2000)
mscolcs <- RunPCA(msclin, npcs = 10)
mscolcs <- FindNeighbors(msclin, dims = 1:10)
mscolcs <- FindClusters(msclin, resolution = 0.2)
mscolcs <- RunUMAP(msclin, dims = 1:10)

bmecs <- ScaleData(bmec)
bmecs <- FindVariableFeatures(bmec, nfeatures = 2000)
bmecs <- RunPCA(bmec, npcs = 7)
bmecs <- FindNeighbors(bmec, dims = 1:7)
bmecs <- FindClusters(bmec, resolution = 0.15)
bmecs <- RunUMAP(bmec, dims = 1:7)
```
<img width="730" alt="gitpic9" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/1ea030ec-5d54-4121-b857-517da0e7848c">

Here you can see the results of our subclustering. Splitting our UMAPs by condition allows us to make preliminary observations of population differences. For example, you can see that in the MSC-lineage cells, there appear to be more cluster 5 cells in the myeloma condition compared to the healthy condition.

Once again, we run a bootstrapping, but this time on the subclusters for the MSC-lineage and BMEC. This supports our first observations and shows these differences are indeed significant. 

![gitpic10](https://github.com/dgagler/5TGM1-Stromal/assets/31450828/db0e430e-12ab-4690-8d2f-b70d8e34ddd1)

Next, we set out to characterize these populations. We want to know what sets these clusters apart, so we will perform a differential gene expression analysis between clusters and visualize the results using a heatmap with the ComplexHeatmap package. The below code shows the process only for the MSC-lineage cells.

```{r}
library(ComplexHeatmap)

# set idents to clusters 
Idents(msclin) <- "seurat_clusters"
# differential expression test (wilcoxon rank-sum by default) between clusters
all.markers <- FindAllMarkers(msclin, assay = "RNA")
all.markers <- all.markers[order(-all.markers$avg_log2FC),] # order by average log-fold change

# we don't want to show all differentially expressed genes in our heatmap so we will extract out the top 15 per cluster
output_gene_order <- c()

# iterate through clusters and get top 15 genes sorted by adjusted p-value and average log-fold change
for(i in c("0", "1", "2", "3", "4", "5")) {
  print(i)

  cluster_genes <- all.markers %>% filter(cluster == i) %>% dplyr::arrange(p_val_adj) %>%
                      dplyr::filter(avg_log2FC > 0) %>%
                      dplyr::select(gene) %>%
                      head(n = 15) %>%
                      purrr::as_vector()

  output_gene_order <- c(output_gene_order, cluster_genes)
}

output_gene_order <- unique(output_gene_order)

# subset our matrix so only include our genes of interest
matrix <- msclin@assays$RNA@scale.data
matrix.roworder.subset <- matrix[rownames(matrix) %in% output_gene_order,]
output_gene_order_sub <- output_gene_order[output_gene_order %in% rownames(matrix.roworder.subset)]

# the ComplexHeatmap function requires you to set your heatmap annotation data beforehand
# we label by experimental condition and seurat cluster and manually set some default R colors
annotation <- HeatmapAnnotation("Disease Status" = msclin@meta.data$condition,
                                "Seurat Cluster" = msclin@meta.data$seurat_clusters,
                                simple_anno_size = unit(3, "mm"),
                                col = list("Disease Status" = c("Healthy" = "green", 
                                                                "Myeloma" = "purple"),
                                           "Seurat Cluster" = c("0" = "#F8766D",
                                                         "1" = "#B79F00",
                                                         "2" = "#00BA38",
                                                         "3" = "#00BFC4",
                                                         "4" = "#619CFF",
                                                         "5" = "#F564E3"
                                                         )))

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# finally, we run the heatmap function
Heatmap(matrix.roworder.subset, cluster_rows = F, cluster_columns = TRUE,
        row_order = output_gene_order_sub,
        show_column_dend = F,
        col = col_fun,
        column_names_gp = gpar(fontsize = 0), row_names_gp = gpar(fontsize = 4.5),
        clustering_distance_columns = "euclidean", clustering_method_columns = "average",
        column_dend_height = unit(40,"mm"), row_dend_width = unit(40, "mm"),
        top_annotation = annotation,
        heatmap_legend_param = list(title = "Scaled Expression")
        )
```
![gitpic11](https://github.com/dgagler/5TGM1-Stromal/assets/31450828/04f447ea-484f-4d93-94e5-2f346c392721)

As mentioned earlier, MSC-lineage cells are known to differentiation into adipocytes and osteo-lineage cells. As such, we will use this framework to characterize our MSC-lineage cells. Marker genes for adipocytes include Adipoq, Lpl, and Mgp. Wif1 has been associated with undifferentiated or multi-potent MSC. Spp1 is a marker gene for early-osteocytes and Bglap is a marker for mature osteocytes, which are bone forming cells. From this heatmap, we can see that cluster 0 is enriched in adipo-genes, cluster 1 expresses adipo genes and Wif1, cluster 2 expresses the early osteocyte gene Spp1, and cluster 3 expresses the mature osteocyte gene Bglap. Cluster 4 expresses similar adipo marker genes to cluster 0 but also expresses a unique set of transcription factors. Cluster 5 is hard to assess from the heatmap, so we generate some violin and feature plots to see if we can glean a new perspective.

![gitpic12](https://github.com/dgagler/5TGM1-Stromal/assets/31450828/d5478752-8fae-4356-a4e9-f3bcb07752d0)

It turns out we can. Interestingly, cluster 5 expresses both adipo and osteo markers and when referring back to our bootstrapping results, we can see that cluster 5 is almost exclusively found in the MM condition and cluster 4 is generally depleted in MM. Taken together, our final subclustering annotation for MSC-lineage cells is as follows:
* cluster 0 = Adipo-MSC
* cluster 1 = Multipotent-MSC
* cluster 2 = Early OLC
* cluster 3 = Mature OLC
* cluster 4 = Adipo-MSC 2
* cluster 5 = MM-OLC

These subclustering results were compiled together with knowledge about MM cell biology to generate a BioRender image which shows an overview of the MSC-lineage populations in our study. Each cluster is shown expressing secreted molecules, adhesion molecules, and receptor molecules.

![gitpic13](https://github.com/dgagler/5TGM1-Stromal/assets/31450828/a6a665a7-f59e-45dc-8237-20b985f65d86)

# Trajectory Analysis

# Gene Set Enrichment Analysis

# EndoMT Analysis

This repository includes the code used for the analysis of this data, including, but not limited to, preprocessing, quality control filtering, integration, clustering, trajectory analysis, and gene set enrichment analysis (GSEA), in addition to the code used to generate both main and supplementary figures.
