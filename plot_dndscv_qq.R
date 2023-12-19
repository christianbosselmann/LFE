## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## This script plots the results of the dNdScv meta-analysis 
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
                 ggpubr,
                 readxl,
                 ggrepel,
                 scales,
                 cowplot)

# helper fns and dict
source("funs.R")
source("dict.R")

# data
df_dndscv <- read_excel_allsheets("data/dndscv.xlsx")
df_dndscv <- rbindlist(df_dndscv, use.names = TRUE, idcol = "phenotype")

### PLOT -----------------------------------------------------------------------
# for the whole model, for each phenotype create a p-value QQ plot and label hits
df_qq <- df_dndscv %>%
  select(phenotype, gene_name, pglobal_cv) 

# fix phenotype annotation
df_qq[df_qq$phenotype == "mcd_α=0.0002304", ]$phenotype <- "MCD"
df_qq[df_qq$phenotype == "leat_α=0.0005376", ]$phenotype <- "LEAT"

# plot parameters
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

# loop over phenotypes of interest
vec_plotphe <- c("MCD", "LEAT")
ls_p <- list()

# plot
for(i in seq_along(vec_plotphe)){
  tmp_phe <- vec_plotphe[[i]]
  
  # censor zero p-values to second-lowest (to avoid 0 issues)
  df_qq[df_qq$pglobal_cv == 0, ]$pglobal_cv <- sort(df_qq[df_qq$pglobal_cv > 0, ]$pglobal_cv, decreasing = FALSE)[1]
  
  tmp <- df_qq %>%
    filter(phenotype == tmp_phe) %>%
    ggplot(aes(sample = -log10(pglobal_cv))) +
    stat_qq(distribution = qexp) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    theme_classic() +
    xlab(log10Pe) +
    ylab(log10Po)
  
  # extract plot data
  df_qq_plot <- ggplot_build(tmp)$data[[1]] %>%
    arrange(desc(y))
  
  # add label - remember to order by pvalue
  df_qq_plot$name <- df_qq %>%
    filter(phenotype == tmp_phe) %>%
    arrange(pglobal_cv) %>%
    pull(gene_name) 
  
  # rebuild
  tmp_p <- df_qq_plot %>%
    mutate(name = ifelse(name %in% vec_known_and_hits, name, NA)) %>% 
    # only label known genes or our candidates
    mutate(name = ifelse(name %in% c(vec_known, "EGFR", "DYRK1A"), name, NA)) %>%
    # color by list
    left_join(df_classification %>% rename(name = Gene)) %>%
    mutate(color = ifelse(Category == "Established lesional epilepsy gene", RColorBrewer::brewer.pal(3, "Set1")[[1]],
                                 ifelse(Category %like% "Established lesional epilepsy gene with first-time statistical support", RColorBrewer::brewer.pal(3, "Set1")[[3]],
                                        ifelse(Category == "Novel lesional epilepsy gene", RColorBrewer::brewer.pal(3, "Set1")[[2]],
                                               "grey")))) %>%
    mutate(color = ifelse(is.na(color), "grey", color)) %>%
    # only label nominally significant genes
    mutate(name = ifelse(y > 2, name, NA))
  
  # plot
  ls_p[[i]] <- tmp_p %>%
    ggplot(aes(theoretical, sample, label = name, color = color)) +
    geom_point(data = tmp_p[tmp_p$color == "grey", ]) +
    geom_point(data = tmp_p[!tmp_p$color == "grey", ]) +
    geom_abline(linetype = "dotted") + 
    geom_label_repel(aes(fontface = "italic"),
                     label.size = .1,
                     min.segment.length = 0, 
                     force = 3,
                     force_pull = .1,
                     label.padding = .1,
                     nudge_x = 0,
                     size = 2,
                     show.legend = FALSE) +
    scale_y_continuous(oob = scales::squish_infinite) +
    scale_x_continuous(oob = scales::squish_infinite) +
    scale_color_identity() +
    theme_classic(base_size = 7) +
    xlab(log10Pe) +
    ylab(log10Po) 
}

p_dndscv_qq <- cowplot::plot_grid(plotlist = ls_p, 
                                  nrow = 2,
                                  align = "hv")
