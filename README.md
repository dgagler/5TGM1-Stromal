# 5TGM1-Stromal Analysis

> This repository contains code used for the analysis of single-cell RNA data as part of [Ghamlouch et al., 2025](https://ehoonline.biomedcentral.com/articles/10.1186/s40164-025-00606-x): *A proinflammatory response and polarized differentiation of stromal elements characterizes the murine myeloma bone marrow niche*.

The raw data for this study can be found on the NCBI's Gene Expression Omnibus (GEO) database under accession [GSE293547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293547). The scripts included in this repository cover the major analytical workflow of the study, outlined as follows:

## Analytical workflow

1. **QC_Filtering.R** - Generation of Seurat objects from single-cell count matrices, demultiplexing of HTOs (and subsequent decision to exclude them due to low quality), QC filtering of cells, and merging of libraries.
2. **Integration_StandardWorkflow.R** - Integration of QC filtered libraries and performing a standard Seurat workflow (PCA, clustering, etc.) on integrated data.
3. **CelltypeAnnotation.R** - Manual annotation of cell populations, identification and reclustering of stromal compartment.
4. **Bootstrapping_Analysis.R** - Bootstrapping analysis to help identify changes in cell type abundance between conditions.
5. **Monocle3_TrajectoryAnalysis.R** - Trajectory analysis via Monocle3.
6. **EndoMT_Analysis.R** - Endothelial-to-mesenchymal (EndoMT) analysis.

The scripts rely on many packages, such as Seurat, ggplot2, and Monocle3, which must be installed prior to utilization of these scripts. Alternatively, these scripts could be imported into something like Rstudio and used interactively.

For questions, please reach out to dylangagler@gmail.com or gareth.morgan@nyulangone.org.
