library(Seurat)
library(Matrix)

cell_subset <- merged

### CONVERTING TO MATRICES/LISTS FOR CELL TYPIST INPUT ###
# Read in Seurat object

# CellTypist requires a log-normalized (to 10,000) expression matrix, a list of cells, and a list of genes/features
raw.matrix <- cell_subset@assays$RNA@data 
cell.names <- Cells(cell_subset@assays$RNA)
gene.names <- rownames(cell_subset@assays$RNA)

# Write out sparse matrix
writeMM(raw.matrix, '/Users/gagled01/morganLab/Waldenstroms/singlecell/WM_PrelimRNA_sparseMatrix.mtx')

# Write out gene and cell lists
write.table(cell.names, '/Users/gagled01/morganLab/Waldenstroms/singlecell/WM_PrelimRNA_cellNames.tsv', quote = F, sep = '\t', 
            col.names = F, row.names = F)
write.table(gene.names, '/Users/gagled01/morganLab/Waldenstroms/singlecell/WM_PrelimRNA_geneNames.tsv', quote = F, sep = '\t',
            col.names = F, row.names = F)
