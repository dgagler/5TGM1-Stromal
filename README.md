# 5TGM1-Stromal Analysis

This repository includes code used for the analysis of single-cell RNA data as part of the Ghamlouch et al., 2025 study in Experimental Hematology & Oncology: *A proinflammatory response and polarized differentiation of stromal elements characterizes the murine myeloma bone marrow niche*.

*Analytical overview is as follows:*

1. **QC_Filtering.R** - Generation of Seurat objects from single-cell count matrices, demultiplexing of HTOs (and subsequent decision to exclude them due to low quality), QC filtering of cells, and merging of libraries.
2. **Integration_StandardWorkflow.R** - Integration of QC filtered libraries and performing a standard Seurat workflow (PCA, clustering, etc.) on integrated data.
3. **CelltypeAnnotation.R** - Manual annotation of cell populations, identification and reclustering of stromal compartment.
4. **Bootstrapping_Analysis.R** - Bootstrapping analysis to help identify changes in cell type abundance between conditions.
5. **Monocle3_TrajectoryAnalysis.R** - Trajectory analysis via Monocle3.
6. **EndoMT_Analysis.R** - Endothelial-to-mesenchymal (EndoMT) analysis.

For questions, please reach out to dylangagler@gmail.com or gareth.morgan@nyulangone.org.
