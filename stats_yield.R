## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## This short script calculates the proportion of cases explained by each gene
##
## Author: Dr. Christian Bosselmann
##
## Date Created: 2023-09-05
##
## Copyright (c) Christian Bosselmann, 2023
## Email: bosselc@ccf.org; christian.bosselmann@gmail.com
##
## -----------------------------------------------------------------------------

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 data.table,
                 readxl)

# data: masterlist
df <- readxl::read_excel("data/MASTERLIST_2023-08-09.xlsx")

### FIRST METHOD ---------------------------------------------------------------
# reduce to columns of interest: just id, phenotype, and somatic genes
df <- df %>%
  select(CCG_ID, category, sub_category, somatic_genes)

# only count cases, not controls
df <- df %>%
  filter(category %in% c("LEAT", "MCD"))

# function: count percentage of samples explained by a gene
explainedByGene <- function(df, gene){
  # if gene is a vector of len > 1, concatenate for grep
  gene <- paste0(gene, collapse = "|")
  nrow(df[grep(gene, df$somatic_genes), ])/nrow(df) 
}

# get list of established LFE genes 
vec_genes <- read_excel("data/known_LFE_genes.xlsx")
vec_genes <- vec_genes %>%
  pull(Gene...6) %>%
  unique
vec_genes <- c(vec_genes, "DYRK1A", "EGFR")

# for each established gene, get percent explained samples
ls_exp <- lapply(vec_genes, explainedByGene, df = df)
names(ls_exp) <- vec_genes

# make pretty
df_exp <- as_tibble(ls_exp) %>%
  pivot_longer(cols = everything(), names_to = "gene") %>%
  arrange(desc(value)) %>%
  mutate(cumsum = cumsum(value))

# plot
p_exp <- df_exp %>%
  mutate(nrow = row_number()) %>%
  mutate(color = ifelse(cumsum < .95*max(df_exp$cumsum), RColorBrewer::brewer.pal(3, "Set1")[[1]], "lightgray")) %>% # color genes that taken together explain 99% of our yield
  ggplot(aes(x = nrow, y = cumsum, label = gene, fill = color)) +
  geom_hline(yintercept = max(df_exp$cumsum)/2, linetype = "dashed", color = "black") +
  geom_col() +
  ggtext::geom_richtext(aes(label = gene), 
                        position = position_stack(vjust = .5),
                        angle = 90, 
                        size = 3,
                        fill = "white", 
                        label.size = 0,
                        fontface = "italic") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(1, 2, 7, 19)) +
  scale_fill_identity() +
  theme_classic() +
  coord_cartesian(expand = FALSE) +
  ylab("Cumulative yield") +
  xlab("Number of genes")

### SECOND METHOD --------------------------------------------------------------
## the first method shown above counts samples once for each gene, which leads to
## a higher cumulative yield than the actual yield; instead, recalculate for 
## gene subsets in order of their respective yield

df_cum <- df_exp %>%
  filter(value > 0) %>%
  select(-cumsum)

# in order of increasing yield, group genes as gene set and get cumulative yield
df_cum$cumsum <- NA
for(i in 1:nrow(df_cum)){
  gene_set <- df_cum[1:i, 1] %>% pull(gene)
  
  # get percentage of samples explained by any gene in the gene set
  df_cum[i, ]$cumsum <- explainedByGene(df, gene_set)
}

# plot
p_cum <- df_cum %>%
  mutate(nrow = row_number()) %>%
  mutate(color = ifelse(cumsum < .95*max(df_cum$cumsum), RColorBrewer::brewer.pal(3, "Set1")[[1]], "lightgray")) %>% # color genes that taken together explain 99% of our yield
  ggplot(aes(x = nrow, y = cumsum, label = gene, fill = color)) +
  geom_hline(yintercept = max(df_cum$cumsum)/2, linetype = "dashed", color = "black") +
  geom_col() +
  ggtext::geom_richtext(aes(label = gene), 
                        position = position_stack(vjust = .5),
                        angle = 90, 
                        size = 2,
                        fill = "white", 
                        label.size = 0,
                        fontface = "italic") +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(pretty(df_cum$cumsum))) +
  scale_x_continuous(breaks = c(1, 2, length(df_cum[df_cum$cumsum < .95*max(df_cum$cumsum),]$cumsum), nrow(df_cum))) +
  scale_fill_identity() +
  theme_classic(base_size = 7) +
  coord_cartesian(expand = FALSE) +
  ylab("Cumulative yield") +
  xlab("Number of genes")

p_cum

# export
pdf("output/cumulative_yield_withnovel.pdf", width = (nrow(df_cum)/12)*3, height = 3)
p_cum
dev.off()
