# Randomly breaking up samples by condition for pseudobulking
# THIS AINT IT

# Case/control breakdown
stromals <- readRDS("/Users/gagled01/morganLab/single-cell/5TG_Mouse_Rerun/Rerun_Mouse5TG_StromalsOnly_AnnotatedObject6_21PCs_NoLib_FIXEDLIBRARIES_6.rds")
healthy <- subset(stromals, subset = case.control == "Healthy")
myeloma <- subset(stromals, subset = case.control == "Myeloma")

# Applying pseudobulk with 3 samples each
healthy@meta.data$pseudobulk <- sample(factor(rep(1:3, length.out = nrow(healthy@meta.data)), labels = paste0("Healthy", 1:3)))
table(healthy@meta.data$pseudobulk)
myeloma@meta.data$pseudobulk <- sample(factor(rep(1:3, length.out = nrow(myeloma@meta.data)), labels = paste0("Myeloma", 1:3)))
table(myeloma@meta.data$pseudobulk)

# Applying pseudobulk with 5 samples each
healthy@meta.data$pseudobulk2 <- sample(factor(rep(1:5, length.out = nrow(healthy@meta.data)), labels = paste0("Healthy", 1:5)))
table(healthy@meta.data$pseudobulk2)
myeloma@meta.data$pseudobulk2 <- sample(factor(rep(1:5, length.out = nrow(myeloma@meta.data)), labels = paste0("Myeloma", 1:5)))
table(myeloma@meta.data$pseudobulk2)

merged.stromals <- merge(healthy, myeloma, merge.data = TRUE) # merge.data = T so we don't wipe out normalized values