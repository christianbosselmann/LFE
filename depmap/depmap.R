## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## Depmap functional analysis
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
                 ggforce,
                 ggrepel)

# helper fns and dict
source("funs.R")
source("dict.R")

### SELECTIVITY/ESSENTIALITY ---------------------------------------------------
# data
df_sel <- readxl::read_excel("depmap/data/fig3-data1-v3.xlsx", sheet = 3)

# plot efficacy vs. selectivity; theta 0.6, X = 1
df_sel <- df_sel %>%
  mutate(label = ifelse(GeneSymbol %in% vec_known_and_hits, GeneSymbol, NA)) %>%
  mutate(color = ifelse(GeneSymbol %in% vec_known_and_hits, "red", "lightgray"))

p_sel <- df_sel %>%
  ggplot(aes(x = q1.efficacy, y = q1.selectivity, label = label, color = color, alpha = color)) +
  geom_point(data = df_sel[df_sel$color == "lightgray", ]) +
  geom_point(data = df_sel[df_sel$color == "red", ]) +
  geom_rug() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = mean(df_sel$q1.efficacy) - 2*sd(df_sel$q1.efficacy), linetype = "dashed") +
  geom_hline(yintercept = mean(df_sel$q1.selectivity) + 2*sd(df_sel$q1.selectivity), linetype = "dashed") +
  geom_label_repel(aes(fontface = "italic"),
                   label.size = .1,
                   min.segment.length = 0,
                   force = 0.1,
                   force_pull = 0.05,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2.5) +
  scale_color_manual(values = c("lightgray", "darkred"), guide = "none") +
  scale_alpha_discrete(range = c(0.2, 1), guide = "none") +
  xlab("Efficacy") +
  ylab("Selectivity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  theme_classic() 

pdf("depmap/sel_eff.pdf", width = 6, height = 4)
p_sel
dev.off()

### DEPENDENCY SCORE DISTRIBUTION ----------------------------------------------
# data
df_dep1 <- readxl::read_excel("depmap/data/fig2-data1-v3.xlsx", sheet = 1) # siRNA
df_dep2 <- readxl::read_excel("depmap/data/fig2-data1-v3.xlsx", sheet = 6) # CRISPR

# preprocess
df_dep1 <- df_dep1 %>%
  mutate_at(vars(contains('ACH')), as.numeric) %>%
  na.omit %>%
  pivot_longer(cols = starts_with("ACH")) %>%
  na.omit %>%
  rename(siRNA = value)

df_dep2 <- df_dep2 %>%
  mutate_at(vars(contains('ACH')), as.numeric) %>%
  na.omit %>%
  pivot_longer(cols = starts_with("ACH")) %>%
  na.omit %>%
  rename(crispr = value)

df_dep <- cbind(df_dep1, df_dep2[,"crispr"])

# plot function
dependencyPlot <- function(df_dep, gene, limit = FALSE){
  
  r_spearman <- cor(df_dep[df_dep$Gene == gene, ]$siRNA,
                    df_dep[df_dep$Gene == gene, ]$crispr,
                    method = "spearman") %>%
    signif(3)
  
  r_pearson <- cor(df_dep[df_dep$Gene == gene, ]$siRNA,
                   df_dep[df_dep$Gene == gene, ]$crispr,
                   method = "pearson") %>%
    signif(3)
  
  dep_label <- paste0("Spearman: ", sprintf('%#.3g', r_spearman), "\n",
                      "Pearson: ", sprintf('%#.3g', r_pearson))
  
  # limit contour plot to known LFE-associated genes
  if(limit){tmp_dep <- df_dep %>% filter(Gene %in% vec_known_and_hits)}else{
    tmp_dep <- df_dep
  }
  
  tmp_dep %>%
    ggplot(aes(x = crispr, y = siRNA)) +
    geom_bin2d(aes(fill = stat(count/max(count))),
               bins = 30, alpha = 1) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(3, "Greys")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(linetype = "dashed") +
    geom_point(data = tmp_dep[tmp_dep$Gene == gene, ]) +
    tune::coord_obs_pred() +
    ggtitle(gene) +
    xlab("CRISPR dependency score") +
    ylab("siRNA dependency score") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(face = "italic")) +
    annotate("text", x = 0.5, y = -Inf, hjust = 0, vjust = -.2, 
             label = dep_label)
}

pls_dep <- lapply(vec_known_and_hits, function(x){dependencyPlot(df_dep, x)})
pls_dep <- lapply(c("DYRK1A", "EGFR", 
                    "BRAF", "MTOR"), function(x){dependencyPlot(df_dep, x)})

# arrange plot
p_dep <- cowplot::plot_grid(plotlist = pls_dep)

# export
pdf("depmap/dep_contour.pdf", height = 8, width = 10)
p_dep
dev.off()

### NETWORK ANALYSIS: FUNCTIONAL -----------------------------------------------
# data: gene cluster membership
df_netw <- readxl::read_excel("depmap/data/fig5-data1-v3.xlsx", sheet = 3)

# plot tSNE subset of known and candidate genes
df_netw <- df_netw %>%
  mutate(is_known = ifelse(GeneSymbol %in% vec_known_and_hits, TRUE, FALSE)) %>%
  mutate(color = ifelse(is_known == TRUE, "Established", "Other")) %>%
  mutate(color = ifelse(GeneSymbol %in% c("DYRK1A", "EGFR"), "Candidate", color))

# contour visualization of established clusters
# this is just done for diagnosis and/or data exploration
p_netw <- df_netw  %>%
  ggplot(aes(x = tSNE_x, y = tSNE_y, color = color, label = GeneSymbol)) +
  geom_point(data = df_netw[df_netw$color == "Other", ]) +
  stat_density_2d(data = df_netw[df_netw$Cluster_Small %in% df_netw[df_netw$color == "Established", ]$Cluster_Small, ], 
                  aes(alpha = ..level.., color = "red"),
                  bins = 5) +
  scale_alpha_continuous(limits = c(0, 5e-4), breaks = 1e-6*seq(0, 4, by = 2))+
  geom_point(data = df_netw[df_netw$color == "Established", ]) +
  geom_point(data = df_netw[df_netw$color == "Candidate", ]) +
  geom_label_repel(data = df_netw[df_netw$GeneSymbol %in% vec_known_and_hits, ],
                   aes(fontface = "italic"),
                   label.size = .1,
                   min.segment.length = 0,
                   force = 0.1,
                   force_pull = 0.05,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 4) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set1")[[2]], 
                                RColorBrewer::brewer.pal(3, "Set1")[[1]], 
                                "gray", "black")) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("tSNE X") +
  ylab("tSNE Y") +
  coord_cartesian(expand = F)

# export
pdf("depmap/dep_netw.pdf", height = 6, width = 6)
p_netw
dev.off()

### NETWORK ANALYSIS: NEIGHBOURHOOD -------------------------------------------- 
## zoom into cluster neighbourhoods of DYRK1A and EGFR for inset plots
# define target window coordinates: DYRK1A
zoom_gene <- "DYRK1A"
zoom_x <- df_netw[df_netw$GeneSymbol == zoom_gene, ]$tSNE_x
zoom_y <- df_netw[df_netw$GeneSymbol == zoom_gene, ]$tSNE_y
zoom_bandwidth <- 25

df_zoom <- df_netw %>%
  mutate(GeneSymbol = ifelse(Cluster_Large == df_netw[df_netw$GeneSymbol == zoom_gene, ]$Cluster_Large, GeneSymbol, NA))

# get gene set for pathway analysis of these neighbourhoods
zoom_genes <- df_zoom %>%
  na.omit %>%
  pull(GeneSymbol)

## enrichR
librarian::shelf(enrichR)

# set up API
websiteLive <- getOption("enrichR.live")
listEnrichrSites()
setEnrichrSite("Enrichr") 

# set list of DBs and select DBs to query
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023",
         "KEGG_2021_Human")

# do query
enriched <- enrichr(zoom_genes, dbs)
enriched <- enriched[dbs]
enriched <- rbindlist(enriched, idcol = "db")

# view
enriched %>%
  filter(Genes %like% zoom_gene) %>%
  filter(Adjusted.P.value < 0.05) %>%
  arrange(Adjusted.P.value) %>%
  view

# labeller function
# note that not all queried ontologies will return >=1 significant association
custom_names <- list(
  "GO_Molecular_Function_2023" = "GO Molecular Function",
  "GO_Biological_Process_2023" = "GO Biological Process",
  "KEGG_2021_Human" = "KEGG"
)

custom_labeller <- function(variable,value){
  return(custom_names[value])
}

# plot
p_enriched <- enriched %>%
  filter(Genes %like% zoom_gene) %>%
  filter(Adjusted.P.value < 0.05) %>%
  ggplot(aes(x = Combined.Score, y = reorder(Term, Combined.Score))) +
  geom_col(position = position_dodge(preserve = "single")) +
  theme_classic(base_size = 7) +
  theme(axis.text.y = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  xlab("Combined Score") +
  ylab("") +
  ggtitle(zoom_gene) +
  ggforce::facet_col(~db, 
             scales = "free",
             space = "free",
             labeller = custom_labeller) +
  coord_cartesian(expand = FALSE)

# repeat everything for EGFR
zoom_gene <- "EGFR"
zoom_x <- df_netw[df_netw$GeneSymbol == zoom_gene, ]$tSNE_x
zoom_y <- df_netw[df_netw$GeneSymbol == zoom_gene, ]$tSNE_y
zoom_bandwidth <- 25

df_zoom <- df_netw %>%
  mutate(GeneSymbol = ifelse(Cluster_Large == df_netw[df_netw$GeneSymbol == zoom_gene, ]$Cluster_Large, GeneSymbol, NA))

# get gene set for pathway analysis of these neighbourhoods
zoom_genes <- df_zoom %>%
  na.omit %>%
  pull(GeneSymbol)

## enrichR
librarian::shelf(enrichR)

# set up API
websiteLive <- getOption("enrichR.live")
listEnrichrSites()
setEnrichrSite("Enrichr") 

# set list of DBs and select DBs to query
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023",
         "KEGG_2021_Human")

# do query
enriched <- enrichr(zoom_genes, dbs)
enriched <- enriched[dbs]
enriched <- rbindlist(enriched, idcol = "db")

# view
enriched %>%
  filter(Genes %like% zoom_gene) %>%
  filter(Adjusted.P.value < 0.05) %>%
  arrange(Adjusted.P.value) %>%
  view

# labeller function
custom_names <- list(
  "GO_Molecular_Function_2023" = "GO Molecular Function",
  "GO_Biological_Process_2023" = "GO Biological Process",
  "KEGG_2021_Human" = "KEGG"
)

custom_labeller <- function(variable,value){
  return(custom_names[value])
}

# plot
p_enriched2 <- enriched %>%
  filter(Genes %like% zoom_gene) %>%
  filter(Adjusted.P.value < 0.05) %>%
  ggplot(aes(x = Combined.Score, y = reorder(Term, Combined.Score))) +
  geom_col(position = position_dodge(preserve = "single")) +
  theme_classic(base_size = 7) +
  theme(axis.text.y = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  xlab("Combined Score") +
  ylab("") +
  ggtitle(zoom_gene) +
  ggforce::facet_col(~db, 
                     scales = "free",
                     space = "free",
                     labeller = custom_labeller) +
  coord_cartesian(expand = FALSE)

# adjust relative plot heights: ratio of terms in plot
rel_height <- c(nrow(p_enriched$data)/(nrow(p_enriched$data)+nrow(p_enriched2$data)), 
  nrow(p_enriched2$data)/(nrow(p_enriched$data)+nrow(p_enriched2$data)))

p_enriched_comb <- cowplot::plot_grid(p_enriched, p_enriched2,
                   rel_heights = rel_height,
                   nrow = 2)

pdf(file = paste0("depmap/enrich_comb", ".pdf"), width = 6, height = 8)
p_enriched_comb
dev.off()
