## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## Depmap functional analysis; pairwise co-dependency
## ref: Gatto et al. 2023; Shimada et al. 2021
##
## Author: Dr. Christian Bosselmann
##
## Date Created: 2023-09-05
##
## Copyright (c) Christian Bosselmann, 2023
## Email: bosselc@ccf.org; christian.bosselmann@gmail.com
##
## -----------------------------------------------------------------------------

### HEADER ---------------------------------------------------------------------
# pkg
library(librarian)
librarian::shelf(tidyverse,
                 data.table,
                 openxlsx,
                 ggExtra,
                 corrplot,
                 ggrepel)

# helper fns and dict
source("funs.R")
source("dict.R")

### DEMETER2 v6 D2-ACHILLES RNAi -----------------------------------------------
# release date 04/20
df_achilles <- read_csv("depmap/data/dependency scores/D2_Achilles_gene_dep_scores.csv") %>%
  rename(gene = `...1`)

# subset to genes of interest
df_achilles <- df_achilles %>%
  separate(gene, into = c("gene", "drop"), sep = "\\s") %>%
  select(-drop)

df_achilles <- df_achilles %>%
  filter(gene %in% vec_known_and_hits)

# median imputation
df_achilles <- df_achilles %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), median(., na.rm = TRUE))))

## pairwise correlation
tmp <- as.matrix(df_achilles[,-1])
rownames(tmp) <- as.vector(df_achilles[,1]$gene)

# flip to genes as columns, cell lines as rows
tmp <- t(tmp)

test_mat <- cor.mtest(tmp, conf.level = 0.95)
cor_mat <- cor(tmp)

# plot and export
pdf("depmap/heatmap_achilles.pdf", height = 6, width = 6)
corrplot(cor_mat,
         method = "color", 
         diag = TRUE,
         type = "lower",
         tl.col = 'black', tl.srt = 45,
         addCoef.col = "darkgray", number.cex = 0.5,
         font = 3,
         number.digits = 1,
         p.mat = test_mat$p,
         sig.level = 0.10,
         insig = "blank",
         cl.ratio = 0.2,
         order = "AOE") # AOE: angular order of eigenvectors
dev.off()

### CRISPR ---------------------------------------------------------------------
# date 05/23
df_crispr <- read_csv("depmap/data/dependency scores/CRISPRGeneEffect.csv")

# transpose
df_crispr <- df_crispr %>%
  pivot_longer(cols = -ModelID, names_to = 'gene') %>% 
  pivot_wider(names_from = ModelID, values_from = value)

# subset to genes of interest
df_crispr <- df_crispr %>%
  separate(gene, into = c("gene", "drop"), sep = "\\s") %>%
  select(-drop)

df_crispr <- df_crispr %>%
  filter(gene %in% vec_known_and_hits)

## pairwise correlation
tmp <- as.matrix(df_crispr[,-1])
rownames(tmp) <- as.vector(df_crispr[,1]$gene)

# flip to genes as columns, cell lines as rows
tmp <- t(tmp)

test_mat <- cor.mtest(tmp, conf.level = 0.95)
cor_mat <- cor(tmp)

# plot and export
pdf("depmap/heatmap_crispr.pdf", height = 6, width = 6)
corrplot(cor_mat,
         method = "color", 
         diag = TRUE,
         type = "lower",
         tl.col = 'black', tl.srt = 45,
         addCoef.col = "darkgray", number.cex = 0.5,
         font = 3,
         number.digits = 1,
         p.mat = test_mat$p,
         sig.level = 0.10,
         insig = "blank",
         cl.ratio = 0.2,
         order = "AOE") # or AOE, angular order of eigenvectors
dev.off()
