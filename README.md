# Single-cell RNA Analysis of the Bone Marrow Stromal Compartment

**Multiple Myeloma (MM)** is a plasma cell malignancy characterized by abnormal antibody production. After receiving driver mutations in the germinal center, the founding myeloma cell travels to the bone marrow where it alters its local cellular environment to be more conducive to its long-term residency and proliferation. This is termed the myeloma niche. Healthy stem cell, which give rise to immune cells in normal contexts, also exist in privileged bone marrow niches. **Stromal cells are an essential part of both healthy stem cell and myeloma cell niches**, providing structural support, growth and survival signals, and regulating immune cell accesse. While it is appreciated that stromal cells support myeloma progression, the specific ways in which the stromal compartment responds to the myeloma cell is not well understood. *This project aims to better understand these changes.*

# Quick Experimental Overview:
<img width="947" alt="image" src="https://github.com/dgagler/5TGM1-Stromal/assets/31450828/ecf3b0ec-f87e-4a62-911c-d3299f522b20">

To study the myeloma niche in a controlled fashion we utilized a mouse model. The 5TGM1 mouse MM model reliably recapitulates both the genetic and physiological features typical of human MM. For this study, 4 mice were injected with 5TGM1 and 4 mice were injected with phosphate buffer saline (PBS) as healthy controls. 1 week after the detection of M-spike proteins in the blood of the 5TGM1 mice the mice were sacrificed and stromal cells were isolated from their bone marrows via flow cytometry and fluoresence-activated cell sorting (FACS). Finally, these cells underwent single-cell RNA sequencing via 10X Genomics. Data was analyzed in R, primarily following a Seurat framework.

# Quality Control
The data gets to me in the form of Cellranger outputs. Cellranger is a software tool which aligns single-cell RNA sequence reads to a reference genome to generate feature-barcode matrices. In the context of single-cell RNA data, features represent genes and barcodes represent individual cells. For this experiment, our cells were divided across 5 libraries, or sequencing runs. We load the libraries as follows:
```{r}
lib1.data <- Read10X('/Users/dgagler/projects/5TGM1_Mouse/sample1/filtered_feature_bc_matrix')
```


# Single-cell RNA Analysis

# Integration

# Cellular Annotation

# Bootstrapping Relative Abundance

# Subcluster Analysis

# Trajectory Analysis

# Gene Set Enrichment Analysis

# EndoMT Analysis

This repository includes the code used for the analysis of this data, including, but not limited to, preprocessing, quality control filtering, integration, clustering, trajectory analysis, and gene set enrichment analysis (GSEA), in addition to the code used to generate both main and supplementary figures.
