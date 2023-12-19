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

### PLOT: DNDS SCATTER PLOT ----------------------------------------------------
df_scatter <- df_dndscv

# fix phenotype annotation
df_scatter[df_scatter$phenotype == "mcd_α=0.0002304", ]$phenotype <- "MCD"
df_scatter[df_scatter$phenotype == "leat_α=0.0005376", ]$phenotype <- "LEAT"

tmp_scatter_mcd <- df_scatter %>%
  filter(phenotype %in% c("MCD")) %>%
  # manual colors
  left_join(df_classification %>% rename(gene_name = Gene)) %>%
  mutate(color = ifelse(Category == "Established lesional epilepsy gene", RColorBrewer::brewer.pal(3, "Set1")[[1]],
                        ifelse(Category %like% "Established lesional epilepsy gene with first-time statistical support", RColorBrewer::brewer.pal(3, "Set1")[[3]],
                               ifelse(Category == "Novel lesional epilepsy gene", RColorBrewer::brewer.pal(3, "Set1")[[2]],
                                      "grey")))) %>%
  mutate(color = ifelse(is.na(color), "grey", color)) %>%
  # sanity check: assert that manually labelled genes are significant for that phenotype
  mutate(color = ifelse(qmis_cv <= 0.05 | qtrunc_cv <= 0.05 | qglobal_cv <= 0.05 | qallsubs_cv <= 0.05, color, "grey")) %>%
  # set label
  mutate(label = ifelse((gene_name %in% vec_known_and_hits) & !color == "grey", gene_name, NA)) 

# plot
p_dnds_scatter_mcd <- tmp_scatter_mcd %>%
  ggplot(aes(x = wmis_cv, y = wnon_cv, color = color, label = label)) +
  geom_point(data = tmp_scatter_mcd[tmp_scatter_mcd$color == "grey", ]) +
  geom_point(data = tmp_scatter_mcd[!tmp_scatter_mcd$color == "grey", ]) +
  geom_label_repel(aes(fontface = "italic"),
                   label.size = .1,
                   min.segment.length = 0, 
                   force = 4,
                   force_pull = .01,
                   label.padding = .1,
                   nudge_x = 0,
                   size = 2,
                   show.legend = FALSE) +
  scale_colour_identity() +
  scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000)) +
  theme_classic(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  xlab("dN/dS missense ratio") +
  ylab("dN/dS nonsense ratio") 

tmp_scatter_leat <- df_scatter %>%
  filter(phenotype %in% c("LEAT")) %>%
  # manual colors
  left_join(df_classification %>% rename(gene_name = Gene)) %>%
  mutate(color = ifelse(Category == "Established lesional epilepsy gene", RColorBrewer::brewer.pal(3, "Set1")[[1]],
                        ifelse(Category %like% "Established lesional epilepsy gene with first-time statistical support", RColorBrewer::brewer.pal(3, "Set1")[[3]],
                               ifelse(Category == "Novel lesional epilepsy gene", RColorBrewer::brewer.pal(3, "Set1")[[2]],
                                      "grey")))) %>%
  mutate(color = ifelse(is.na(color), "grey", color)) %>%
  # sanity check: assert that manually labelled genes are significant for that phenotype
  mutate(color = ifelse(qmis_cv <= 0.05 | qtrunc_cv <= 0.05 | qglobal_cv <= 0.05 | qallsubs_cv <= 0.05, color, "grey")) %>%
  # set label
  mutate(label = ifelse((gene_name %in% vec_known_and_hits) & !color == "grey", gene_name, NA)) 

# plot
p_dnds_scatter_leat <- tmp_scatter_leat %>%
  ggplot(aes(x = wmis_cv, y = wnon_cv, color = color, label = label)) +
  geom_point(data = tmp_scatter_leat[tmp_scatter_leat$color == "grey", ]) +
  geom_point(data = tmp_scatter_leat[!tmp_scatter_leat$color == "grey", ]) +
  geom_label_repel(aes(fontface = "italic"),
                   label.size = .1,
                   min.segment.length = 0, 
                   force = 3,
                   force_pull = .1,
                   label.padding = .1,
                   nudge_x = 0,
                   size = 2,
                   show.legend = FALSE) +
  scale_colour_identity() +
  scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000)) +
  theme_classic(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  xlab("dN/dS missense ratio") +
  ylab("dN/dS nonsense ratio") 

p_dndscv_scatter <- cowplot::plot_grid(p_dnds_scatter_mcd, 
                                       p_dnds_scatter_leat, 
                                       nrow = 2,
                                       align = "hv")

# export via collect_plots.R