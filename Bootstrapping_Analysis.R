rm(list = ls()) # clear environment
options(warn=-1) # turn off warning messages globally

library(Seurat)
library(dplyr)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(patchwork)

options(future.globals.maxSize=1000000000000000) # Set max global size so we don't run out of memory

stromals <- readRDS("/Users/gagled01/morganLab/single-cell/5TGM1_Stromal/objects/Rerun_Mouse5TG_StromalsOnly_AnnotatedObject6_21PCs_NoLib_FIXEDLIBRARIES_6.rds")
myeloma <- subset(stromals, subset = case.control %in% "Myeloma")
healthy <- subset(stromals, subset = case.control %in% "Healthy")

# Total cell counts by condition
myeloma.abundances <- myeloma@meta.data$celltypes_ECsubbed
healthy.abundances <- healthy@meta.data$celltypes_ECsubbed

# Relative abundances by condition
#null.relabuns <- table(stromals@meta.data$celltypes_ECsubbed)/sum(table(stromals@meta.data$celltypes_ECsubbed))
myeloma.relabuns <- table(myeloma@meta.data$celltypes_ECsubbed)/sum(table(myeloma@meta.data$celltypes_ECsubbed))
healthy.relabuns <- table(healthy@meta.data$celltypes_ECsubbed)/sum(table(healthy@meta.data$celltypes_ECsubbed))

# Isolate values for each observed cell type relative abundance for MM
mm.observed.sec <- myeloma.relabuns[1]
mm.observed.fibro <- myeloma.relabuns[2]
mm.observed.msc <- myeloma.relabuns[3]
mm.observed.olc <- myeloma.relabuns[4]
mm.observed.peri <- myeloma.relabuns[5]
mm.observed.aec <- myeloma.relabuns[6]
mm.observed.chondro <- myeloma.relabuns[7]

# Isolate values for each observed cell type relative abundance for Healthy
healthy.observed.sec <- healthy.relabuns[1]
healthy.observed.fibro <- healthy.relabuns[2]
healthy.observed.msc <- healthy.relabuns[3]
healthy.observed.olc <- healthy.relabuns[4]
healthy.observed.peri <- healthy.relabuns[5]
healthy.observed.aec <- healthy.relabuns[6]
healthy.observed.chondro <- healthy.relabuns[7]

# Total cells by condition
n_cells_myeloma <- ncol(myeloma)
n_cells_healthy <- ncol(healthy)


# Randomly subsample with replacement and spit out relative abundances
myeloma.resamples <- lapply(1:100, function(i) table(sample(myeloma.abundances, replace = T, prob = ))/n_cells_myeloma)
healthy.resamples <- lapply(1:100, function(i) table(sample(healthy.abundances, replace = T, prob = ))/n_cells_healthy)

# Relative abundance distributions for each celltype -- MM
mm.resamples.sec <- sapply(myeloma.resamples,"[[",1)
mm.resamples.fibro <- sapply(myeloma.resamples,"[[",2)
mm.resamples.msc <- sapply(myeloma.resamples,"[[",3)
mm.resamples.olc <- sapply(myeloma.resamples,"[[",4)
mm.resamples.peri <- sapply(myeloma.resamples,"[[",5)
mm.resamples.aec <- sapply(myeloma.resamples,"[[",6)
mm.resamples.chondro <- sapply(myeloma.resamples,"[[",7)

# Relative abundance distributions for each celltype -- Healthy
healthy.resamples.sec <- sapply(healthy.resamples,"[[",1)
healthy.resamples.fibro <- sapply(healthy.resamples,"[[",2)
healthy.resamples.msc <- sapply(healthy.resamples,"[[",3)
healthy.resamples.olc <- sapply(healthy.resamples,"[[",4)
healthy.resamples.peri <- sapply(healthy.resamples,"[[",5)
healthy.resamples.aec <- sapply(healthy.resamples,"[[",6)
healthy.resamples.chondro <- sapply(healthy.resamples,"[[",7)

sec.plotting.df <- melt(data.frame(Myeloma = mm.resamples.sec, Healthy = healthy.resamples.sec, Celltype = "SEC"))
colnames(sec.plotting.df)[2] <- "Condition"
fibro.plotting.df <- melt(data.frame(Myeloma = mm.resamples.fibro, Healthy = healthy.resamples.fibro, Celltype = "Fibroblasts"))
colnames(fibro.plotting.df)[2] <- "Condition"
msc.plotting.df <- melt(data.frame(Myeloma = mm.resamples.msc, Healthy = healthy.resamples.msc, Celltype = "MSC"))
colnames(msc.plotting.df)[2] <- "Condition"
olc.plotting.df <- melt(data.frame(Myeloma = mm.resamples.olc, Healthy = healthy.resamples.olc, Celltype = "OLC"))
colnames(olc.plotting.df)[2] <- "Condition"
peri.plotting.df <- melt(data.frame(Myeloma = mm.resamples.peri, Healthy = healthy.resamples.peri, Celltype = "Pericytes"))
colnames(peri.plotting.df)[2] <- "Condition"
aec.plotting.df <- melt(data.frame(Myeloma = mm.resamples.aec, Healthy = healthy.resamples.aec, Celltype = "AEC"))
colnames(aec.plotting.df)[2] <- "Condition"
chondro.plotting.df <- melt(data.frame(Myeloma = mm.resamples.chondro, Healthy = healthy.resamples.chondro, Celltype = "Chondrocytes"))
colnames(chondro.plotting.df)[2] <- "Condition"

full.plotting.df <- do.call("rbind", list(sec.plotting.df, fibro.plotting.df, msc.plotting.df, olc.plotting.df, peri.plotting.df, aec.plotting.df, chondro.plotting.df))

library(rstatix)
sec.stat.test <- sec.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
fibro.stat.test <- fibro.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
msc.stat.test <- msc.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
olc.stat.test <- olc.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
peri.stat.test <- peri.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
aec.stat.test <- aec.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
chondro.stat.test <- chondro.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()

msc.other.ttest <- t.test(msc.plotting.df$value ~ as.character(msc.plotting.df$Condition), alternative = c("two.sided"))
olc.other.ttest <- t.test(olc.plotting.df$value ~ as.character(olc.plotting.df$Condition), alternative = c("two.sided"))
sec.other.ttest <- t.test(sec.plotting.df$value ~ as.character(sec.plotting.df$Condition), alternative = c("two.sided"))
aec.other.ttest <- t.test(aec.plotting.df$value ~ as.character(aec.plotting.df$Condition), alternative = c("two.sided"))
fibro.other.ttest <- t.test(fibro.plotting.df$value ~ as.character(fibro.plotting.df$Condition), alternative = c("two.sided"))
chondro.other.ttest <- t.test(chondro.plotting.df$value ~ as.character(chondro.plotting.df$Condition), alternative = c("two.sided"))
peri.other.ttest <- t.test(peri.plotting.df$value ~ as.character(peri.plotting.df$Condition), alternative = c("two.sided"))

msc.lower.conf <- msc.other.ttest$conf.int[1]
msc.upper.conf <- msc.other.ttest$conf.int[2]

olc.lower.conf <- olc.other.ttest$conf.int[1]
olc.upper.conf <- olc.other.ttest$conf.int[2]

sec.lower.conf <- sec.other.ttest$conf.int[1]
sec.upper.conf <- sec.other.ttest$conf.int[2]

aec.lower.conf <- aec.other.ttest$conf.int[1]
aec.upper.conf <- aec.other.ttest$conf.int[2]

fibro.lower.conf <- fibro.other.ttest$conf.int[1]
fibro.upper.conf <- fibro.other.ttest$conf.int[2]

chondro.lower.conf <- chondro.other.ttest$conf.int[1]
chondro.upper.conf <- chondro.other.ttest$conf.int[2]

peri.lower.conf <- peri.other.ttest$conf.int[1]
peri.upper.conf <- peri.other.ttest$conf.int[2]


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
#############################################################
fibro.bxp <- ggboxplot(fibro.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE) +
  ylab("") + xlab("") + ggtitle("Fibroblasts") +
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + NoLegend()
fibro.stat.test <- fibro.stat.test %>% add_xy_position("Condition")
fibro.bxp <- fibro.bxp + 
  stat_pvalue_manual(fibro.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
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
#############################################################
olc.bxp <- ggboxplot(olc.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE) +
  ylab("") + xlab("") + ggtitle("OLC") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + NoLegend()
olc.stat.test <- olc.stat.test %>% add_xy_position("Condition")
olc.bxp <- olc.bxp + 
  stat_pvalue_manual(olc.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
peri.bxp <- ggboxplot(peri.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE) +
  ylab("") + xlab("") + ggtitle("Pericytes") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) # Keeping legend for this one
peri.stat.test <- peri.stat.test %>% add_xy_position("Condition")
peri.bxp <- peri.bxp + 
  stat_pvalue_manual(peri.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
aec.bxp <- ggboxplot(aec.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE) +
  ylab("") + xlab("") + ggtitle("AEC") + 
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + NoLegend()
aec.stat.test <- aec.stat.test %>% add_xy_position("Condition")
aec.bxp <- aec.bxp + 
  stat_pvalue_manual(aec.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
chondro.bxp <- ggboxplot(chondro.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE) +
  ylab("") + xlab("") + ggtitle("Chondrocytes") + 
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + NoLegend()
chondro.stat.test <- chondro.stat.test %>% add_xy_position("Condition")
chondro.bxp <- chondro.bxp + 
  stat_pvalue_manual(chondro.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

msc.bxp | olc.bxp | sec.bxp | aec.bxp | fibro.bxp | chondro.bxp | peri.bxp
ggsave("/Users/gagled01/morganLab/single-cell/5TGM1_Stromal/figures/100x_Bootstrap_CelltypeBoxplots.png", height = 5, width = 20)

#bmecs <- readRDS("/Users/gagled01/morganLab/single-cell/5TGM1_Stromal/objects/5TGM1_BMEC_Object.rds")
bmecs.myeloma <- subset(bmecs, subset = condition %in% "5TGM1")
bmecs.healthy <- subset(bmecs, subset = condition %in% "PBS")

# Total cell counts by condition
bmec.myeloma.abundances <- bmecs.myeloma@meta.data$seurat_clusters
bmec.healthy.abundances <- bmecs.healthy@meta.data$seurat_clusters

# Relative abundances by condition
bmec.myeloma.relabuns <- table(bmecs.myeloma@meta.data$seurat_clusters)/sum(table(bmecs.myeloma@meta.data$seurat_clusters))
bmec.healthy.relabuns <- table(bmecs.healthy@meta.data$seurat_clusters)/sum(table(bmecs.healthy@meta.data$seurat_clusters))

# Isolate values for each observed cell type relative abundance for MM
bmec.mm.observed.c0 <-  bmec.myeloma.relabuns[1]
bmec.mm.observed.c1 <-  bmec.myeloma.relabuns[2]
bmec.mm.observed.c2 <-  bmec.myeloma.relabuns[3]
bmec.mm.observed.c3 <-  bmec.myeloma.relabuns[4]
bmec.mm.observed.c4 <-  bmec.myeloma.relabuns[5]
bmec.mm.observed.c5 <-  bmec.myeloma.relabuns[6]

# Isolate values for each observed cell type relative abundance for Healthy
bmec.healthy.observed.c0 <- bmec.healthy.relabuns[1]
bmec.healthy.observed.c1 <- bmec.healthy.relabuns[2]
bmec.healthy.observed.c2 <- bmec.healthy.relabuns[3]
bmec.healthy.observed.c3 <- bmec.healthy.relabuns[4]
bmec.healthy.observed.c4 <- bmec.healthy.relabuns[5]
bmec.healthy.observed.c5 <- bmec.healthy.relabuns[6]

# Total cells by condition
bmec.n_cells_myeloma <- ncol(bmecs.myeloma)
bmec.n_cells_healthy <- ncol(bmecs.healthy)

# Randomly subsample with replacement and spit out relative abundances
bmec.myeloma.resamples <- lapply(1:100, function(i) table(sample(bmec.myeloma.abundances, replace = T))/bmec.n_cells_myeloma)
bmec.healthy.resamples <- lapply(1:100, function(i) table(sample(bmec.healthy.abundances, replace = T))/bmec.n_cells_healthy)

# Relative abundance distributions for each celltype -- MM
bmec.mm.resamples.c0 <- sapply(bmec.myeloma.resamples,"[[",1)
bmec.mm.resamples.c1 <- sapply(bmec.myeloma.resamples,"[[",2)
bmec.mm.resamples.c2 <- sapply(bmec.myeloma.resamples,"[[",3)
bmec.mm.resamples.c3 <- sapply(bmec.myeloma.resamples,"[[",4)
bmec.mm.resamples.c4 <- sapply(bmec.myeloma.resamples,"[[",5)
bmec.mm.resamples.c5 <- sapply(bmec.myeloma.resamples,"[[",6)

# Relative abundance distributions for each celltype -- Healthy
bmec.healthy.resamples.c0 <- sapply(bmec.healthy.resamples,"[[",1)
bmec.healthy.resamples.c1 <- sapply(bmec.healthy.resamples,"[[",2)
bmec.healthy.resamples.c2 <- sapply(bmec.healthy.resamples,"[[",3)
bmec.healthy.resamples.c3 <- sapply(bmec.healthy.resamples,"[[",4)
bmec.healthy.resamples.c4 <- sapply(bmec.healthy.resamples,"[[",5)
bmec.healthy.resamples.c5 <- sapply(bmec.healthy.resamples,"[[",6)

bmec.c0.plotting.df <- melt(data.frame(Myeloma = bmec.mm.resamples.c0, Healthy = bmec.healthy.resamples.c0, Cluster = "0"))
colnames(bmec.c0.plotting.df)[2] <- "Condition"
bmec.c1.plotting.df <- melt(data.frame(Myeloma = bmec.mm.resamples.c1, Healthy = bmec.healthy.resamples.c1, Cluster = "1"))
colnames(bmec.c1.plotting.df)[2] <- "Condition"
bmec.c2.plotting.df <- melt(data.frame(Myeloma = bmec.mm.resamples.c2, Healthy = bmec.healthy.resamples.c2, Cluster = "2"))
colnames(bmec.c2.plotting.df)[2] <- "Condition"
bmec.c3.plotting.df <- melt(data.frame(Myeloma = bmec.mm.resamples.c3, Healthy = bmec.healthy.resamples.c3, Cluster = "3"))
colnames(bmec.c3.plotting.df)[2] <- "Condition"
bmec.c4.plotting.df <- melt(data.frame(Myeloma = bmec.mm.resamples.c4, Healthy = bmec.healthy.resamples.c4, Cluster = "4"))
colnames(bmec.c4.plotting.df)[2] <- "Condition"
bmec.c5.plotting.df <- melt(data.frame(Myeloma = bmec.mm.resamples.c5, Healthy = bmec.healthy.resamples.c5, Cluster = "5"))
colnames(bmec.c5.plotting.df)[2] <- "Condition"

bmec.c0.stat.test <- bmec.c0.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
bmec.c1.stat.test <- bmec.c1.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
bmec.c2.stat.test <- bmec.c2.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
bmec.c3.stat.test <- bmec.c3.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
bmec.c4.stat.test <- bmec.c4.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
bmec.c5.stat.test <- bmec.c5.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()

# Reorder factor
bmec.c0.plotting.df$Condition <- factor(bmec.c0.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c1.plotting.df$Condition <- factor(bmec.c1.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c2.plotting.df$Condition <- factor(bmec.c2.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c3.plotting.df$Condition <- factor(bmec.c3.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c4.plotting.df$Condition <- factor(bmec.c4.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c5.plotting.df$Condition <- factor(bmec.c5.plotting.df$Condition, levels = c("Healthy", "Myeloma"))

###
bmec.c0.bxp <- ggboxplot(bmec.c0.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("Relative Abundance") + xlab("") + ggtitle("Cluster 0") +
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
bmec.c0.stat.test <- bmec.c0.stat.test %>% add_xy_position("Condition")
bmec.c0.bxp <- bmec.c0.bxp + 
  stat_pvalue_manual(bmec.c0.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
bmec.c1.plotting.df$Condition <- factor(bmec.c1.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c1.bxp <- ggboxplot(bmec.c1.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("") + xlab("") + ggtitle("Cluster 1") +
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
bmec.c1.stat.test <- bmec.c1.stat.test %>% add_xy_position("Condition")
bmec.c1.bxp <- bmec.c1.bxp + 
  stat_pvalue_manual(bmec.c1.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
bmec.c2.plotting.df$Condition <- factor(bmec.c2.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c2.bxp <- ggboxplot(bmec.c2.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("") + xlab("") + ggtitle("Cluster 2") +
  theme_bw(base_size = 16) +
theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
bmec.c2.stat.test <- bmec.c2.stat.test %>% add_xy_position("Condition")
bmec.c2.bxp <- bmec.c2.bxp + 
  stat_pvalue_manual(bmec.c2.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
bmec.c3.plotting.df$Condition <- factor(bmec.c3.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c3.bxp <- ggboxplot(bmec.c3.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("Relative Abundance") + xlab("") + ggtitle("Cluster 3") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
bmec.c3.stat.test <- bmec.c3.stat.test %>% add_xy_position("Condition")
bmec.c3.bxp <- bmec.c3.bxp + 
  stat_pvalue_manual(bmec.c3.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
bmec.c4.plotting.df$Condition <- factor(bmec.c4.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c4.bxp <- ggboxplot(bmec.c4.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("") + xlab("") + ggtitle("Cluster 4") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
bmec.c4.stat.test <- bmec.c4.stat.test %>% add_xy_position("Condition")
bmec.c4.bxp <- bmec.c4.bxp + 
  stat_pvalue_manual(bmec.c4.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
bmec.c5.plotting.df$Condition <- factor(bmec.c5.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
bmec.c5.bxp <- ggboxplot(bmec.c5.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("") + xlab("") + ggtitle("Cluster 5") + 
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
bmec.c5.stat.test <- bmec.c5.stat.test %>% add_xy_position("Condition")
bmec.c5.bxp <- bmec.c5.bxp + 
  stat_pvalue_manual(bmec.c5.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

(bmec.c0.bxp | bmec.c1.bxp | bmec.c2.bxp) / (bmec.c3.bxp | bmec.c4.bxp | bmec.c5.bxp)
ggsave("/Users/gagled01/morganLab/5TGM1_Stromal/figures/100x_Bootstrap_BMEC_Cluster_Boxplots.png", height = 7, width = 8)

#msclin <- readRDS("/Users/gagled01/morganLab/single-cell/5TGM1_Stromal/objects/5TGM1_MSCOLC_Object.rds")
msclins.myeloma <- subset(msclin, subset = condition %in% "5TGM1")
msclins.healthy <- subset(msclin, subset = condition %in% "PBS")

# Total cell counts by condition
msclin.myeloma.abundances <- msclins.myeloma@meta.data$seurat_clusters
msclin.healthy.abundances <- msclins.healthy@meta.data$seurat_clusters

# Relative abundances by condition
msclin.myeloma.relabuns <- table(msclins.myeloma@meta.data$seurat_clusters)/sum(table(msclins.myeloma@meta.data$seurat_clusters))
msclin.healthy.relabuns <- table(msclins.healthy@meta.data$seurat_clusters)/sum(table(msclins.healthy@meta.data$seurat_clusters))

# Isolate values for each observed cell type relative abundance for MM
msclin.mm.observed.c0 <-  msclin.myeloma.relabuns[1]
msclin.mm.observed.c1 <-  msclin.myeloma.relabuns[2]
msclin.mm.observed.c2 <-  msclin.myeloma.relabuns[3]
msclin.mm.observed.c3 <-  msclin.myeloma.relabuns[4]
msclin.mm.observed.c4 <-  msclin.myeloma.relabuns[5]
msclin.mm.observed.c5 <-  msclin.myeloma.relabuns[6]

# Isolate values for each observed cell type relative abundance for Healthy
msclin.healthy.observed.c0 <- msclin.healthy.relabuns[1]
msclin.healthy.observed.c1 <- msclin.healthy.relabuns[2]
msclin.healthy.observed.c2 <- msclin.healthy.relabuns[3]
msclin.healthy.observed.c3 <- msclin.healthy.relabuns[4]
msclin.healthy.observed.c4 <- msclin.healthy.relabuns[5]
msclin.healthy.observed.c5 <- msclin.healthy.relabuns[6]

# Total cells by condition
msclin.n_cells_myeloma <- ncol(msclins.myeloma)
msclin.n_cells_healthy <- ncol(msclins.healthy)

# Randomly subsample with replacement and spit out relative abundances
msclin.myeloma.resamples <- lapply(1:100, function(i) table(sample(msclin.myeloma.abundances, replace = T))/msclin.n_cells_myeloma)
msclin.healthy.resamples <- lapply(1:100, function(i) table(sample(msclin.healthy.abundances, replace = T))/msclin.n_cells_healthy)

# Relative abundance distributions for each celltype -- MM
msclin.mm.resamples.c0 <- sapply(msclin.myeloma.resamples,"[[",1)
msclin.mm.resamples.c1 <- sapply(msclin.myeloma.resamples,"[[",2)
msclin.mm.resamples.c2 <- sapply(msclin.myeloma.resamples,"[[",3)
msclin.mm.resamples.c3 <- sapply(msclin.myeloma.resamples,"[[",4)
msclin.mm.resamples.c4 <- sapply(msclin.myeloma.resamples,"[[",5)
msclin.mm.resamples.c5 <- sapply(msclin.myeloma.resamples,"[[",6)

# Relative abundance distributions for each celltype -- Healthy
msclin.healthy.resamples.c0 <- sapply(msclin.healthy.resamples,"[[",1)
msclin.healthy.resamples.c1 <- sapply(msclin.healthy.resamples,"[[",2)
msclin.healthy.resamples.c2 <- sapply(msclin.healthy.resamples,"[[",3)
msclin.healthy.resamples.c3 <- sapply(msclin.healthy.resamples,"[[",4)
msclin.healthy.resamples.c4 <- sapply(msclin.healthy.resamples,"[[",5)
msclin.healthy.resamples.c5 <- sapply(msclin.healthy.resamples,"[[",6)

msclin.c0.plotting.df <- melt(data.frame(Myeloma = msclin.mm.resamples.c0, Healthy = msclin.healthy.resamples.c0, Cluster = "0"))
colnames(msclin.c0.plotting.df)[2] <- "Condition"
msclin.c1.plotting.df <- melt(data.frame(Myeloma = msclin.mm.resamples.c1, Healthy = msclin.healthy.resamples.c1, Cluster = "1"))
colnames(msclin.c1.plotting.df)[2] <- "Condition"
msclin.c2.plotting.df <- melt(data.frame(Myeloma = msclin.mm.resamples.c2, Healthy = msclin.healthy.resamples.c2, Cluster = "2"))
colnames(msclin.c2.plotting.df)[2] <- "Condition"
msclin.c3.plotting.df <- melt(data.frame(Myeloma = msclin.mm.resamples.c3, Healthy = msclin.healthy.resamples.c3, Cluster = "3"))
colnames(msclin.c3.plotting.df)[2] <- "Condition"
msclin.c4.plotting.df <- melt(data.frame(Myeloma = msclin.mm.resamples.c4, Healthy = msclin.healthy.resamples.c4, Cluster = "4"))
colnames(msclin.c4.plotting.df)[2] <- "Condition"
msclin.c5.plotting.df <- melt(data.frame(Myeloma = msclin.mm.resamples.c5, Healthy = msclin.healthy.resamples.c5, Cluster = "5"))
colnames(msclin.c5.plotting.df)[2] <- "Condition"

msclin.c0.stat.test <- msclin.c0.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
msclin.c1.stat.test <- msclin.c1.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
msclin.c2.stat.test <- msclin.c2.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
msclin.c3.stat.test <- msclin.c3.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
msclin.c4.stat.test <- msclin.c4.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()
msclin.c5.stat.test <- msclin.c5.plotting.df %>%
  t_test(value ~ Condition) %>%
  add_significance()

# Reorder factors
msclin.c0.plotting.df$Condition <- factor(msclin.c0.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
msclin.c1.plotting.df$Condition <- factor(msclin.c1.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
msclin.c2.plotting.df$Condition <- factor(msclin.c2.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
msclin.c3.plotting.df$Condition <- factor(msclin.c3.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
msclin.c4.plotting.df$Condition <- factor(msclin.c4.plotting.df$Condition, levels = c("Healthy", "Myeloma"))
msclin.c5.plotting.df$Condition <- factor(msclin.c5.plotting.df$Condition, levels = c("Healthy", "Myeloma"))

msclin.c0.bxp <- ggboxplot(msclin.c0.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("Relative Abundance") + xlab("") + ggtitle("Cluster 0") +
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
msclin.c0.stat.test <- msclin.c0.stat.test %>% add_xy_position("Condition")
msclin.c0.bxp <- msclin.c0.bxp + 
  stat_pvalue_manual(msclin.c0.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
msclin.c1.bxp <- ggboxplot(msclin.c1.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("") + xlab("") + ggtitle("Cluster 1") +
  theme_bw(base_size = 16) + 
  theme(
        axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
msclin.c1.stat.test <- msclin.c1.stat.test %>% add_xy_position("Condition")
msclin.c1.bxp <- msclin.c1.bxp + 
  stat_pvalue_manual(msclin.c1.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
msclin.c2.bxp <- ggboxplot(msclin.c2.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("") + xlab("") + ggtitle("Cluster 2") +
  theme_bw(base_size = 16) +
theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
msclin.c2.stat.test <- msclin.c2.stat.test %>% add_xy_position("Condition")
msclin.c2.bxp <- msclin.c2.bxp + 
  stat_pvalue_manual(msclin.c2.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
msclin.c3.bxp <- ggboxplot(msclin.c3.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("Relative Abundance") + xlab("") + ggtitle("Cluster 3") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
msclin.c3.stat.test <- msclin.c3.stat.test %>% add_xy_position("Condition")
msclin.c3.bxp <- msclin.c3.bxp + 
  stat_pvalue_manual(msclin.c3.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
msclin.c4.bxp <- ggboxplot(msclin.c4.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("") + xlab("") + ggtitle("Cluster 4") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
msclin.c4.stat.test <- msclin.c4.stat.test %>% add_xy_position("Condition")
msclin.c4.bxp <- msclin.c4.bxp + 
  stat_pvalue_manual(msclin.c4.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#############################################################
msclin.c5.bxp <- ggboxplot(msclin.c5.plotting.df, x = "Condition", y = "value", fill = "Condition", notch = TRUE, add = "jitter", add.params = list(size = 1, alpha = 0.75)) +
  ylab("") + xlab("") + ggtitle("Cluster 5") + 
  theme_bw(base_size = 16) + 
  theme(
        axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face="bold")) + NoLegend()
msclin.c5.stat.test <- msclin.c5.stat.test %>% add_xy_position("Condition")
msclin.c5.bxp <- msclin.c5.bxp + 
  stat_pvalue_manual(msclin.c5.stat.test, label = "T-test, p = {p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

(msclin.c0.bxp | msclin.c1.bxp | msclin.c2.bxp) / (msclin.c3.bxp | msclin.c4.bxp | msclin.c5.bxp)
ggsave("/Users/gagled01/morganLab/5TGM1_Stromal/figures/100x_Bootstrap_MSClin_Cluster_Boxplots.png", height = 7, width = 8)
