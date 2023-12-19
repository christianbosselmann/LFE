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
# max. Q value to plot
label_threshold <- 0.05

# labeller
plot_names <- c(
  "all_α" = "All samples",
  "mcd_α" = "MCD",
  "fcd_i_α" = "FCD I",
  "fcd_ii_α" = "FCD II",
  "moghe_α" = "MOGHE",
  "leat_α" = "LEAT",
  "gg_α" = "GG",
  "dnet_α" = "DNET",
  "other_leat_α" = "Other LEAT"
)

# preprocess dataframe
tmp <- df_dndscv %>%
  separate(phenotype, c("phenotype", "threshold"), "=", extra = "merge") %>%
  mutate(threshold = as.numeric(threshold))

# manual colors
tmp <- tmp %>%
  left_join(df_classification %>% rename(gene_name = Gene)) %>%
  mutate(significance = ifelse(Category == "Established lesional epilepsy gene", RColorBrewer::brewer.pal(3, "Set1")[[1]],
                               ifelse(Category %like% "Established lesional epilepsy gene with first-time statistical support", RColorBrewer::brewer.pal(3, "Set1")[[3]],
                                      ifelse(Category == "Novel lesional epilepsy gene", RColorBrewer::brewer.pal(3, "Set1")[[2]],
                                             "grey")))) %>%
  mutate(significance = ifelse(is.na(significance), "grey", significance))


p_dndscv <- tmp %>%
  mutate(phenotype = factor(phenotype, levels = names(plot_names))) %>%
  ggplot(aes(x = or_case_ctrl, y = qglobal_cv, color = significance)) +
  geom_point() +
  geom_hline(yintercept = label_threshold, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_label_repel(data = . %>% group_by(phenotype) %>% slice_min(qglobal_cv, n = 99, with_ties = F) %>% filter(or_case_ctrl > 0 & qglobal_cv < label_threshold), 
                   aes(label = sprintf('%s', gene_name), fontface = "italic"), 
                   label.size = .1,
                   min.segment.length = 0, 
                   label.padding = .1,
                   force = 1,
                   force_pull = .5,
                   size = 2,
                   show.legend = FALSE) +
  scale_y_continuous(trans = reverselog_trans(10), 
                     breaks = c(1e0, label_threshold, 1e-5, 1e-10, 1e-15),
                     oob = scales::oob_squish_infinite,
                     labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(trans = "log10", 
                     oob = scales::oob_squish_infinite,
                     labels = function(x) ifelse(x == 0, "0", x)) +
  scale_colour_identity() +
  theme_classic(base_size = 7) +
  xlab("OR Case vs. Control") +
  ylab("dNdScv Q-value") +
  coord_cartesian(xlim = c(-10, 100), ylim = c(1, 1e-15), clip = "on") +
  facet_wrap(phenotype ~ ., labeller = as_labeller(plot_names), nrow = 2) +
  theme(strip.background = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA))

### PLOT: SUBSET ---------------------------------------------------------------
# Just data from MCD and LEAT samples, respectively
p_dndscv_mcd <- tmp %>%
  filter(phenotype == "mcd_α") %>%
  ggplot(aes(x = or_case_ctrl, y = qglobal_cv, color = significance)) +
  geom_point() +
  geom_hline(yintercept = label_threshold, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_label_repel(data = . %>% group_by(phenotype) %>% slice_min(qglobal_cv, n = 99, with_ties = F) %>% filter(or_case_ctrl > 0 & qglobal_cv < label_threshold), 
                   aes(label = sprintf('%s', gene_name), fontface = "italic"), 
                   label.size = .1,
                   min.segment.length = 0, 
                   label.padding = .1,
                   force = 1,
                   force_pull = .5,
                   size = 2,
                   show.legend = FALSE) +
  scale_y_continuous(trans = reverselog_trans(10), 
                     breaks = c(1e0, label_threshold, 1e-5, 1e-10, 1e-15),
                     oob = scales::oob_squish_infinite,
                     labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(trans = "log10", 
                     oob = scales::oob_squish_infinite,
                     labels = function(x) ifelse(x == 0, "0", x)) +
  scale_colour_identity() +
  theme_classic(base_size = 7) +
  xlab("OR Case vs. Control") +
  ylab("dNdScv Q-value") +
  coord_cartesian(xlim = c(-10, 100), ylim = c(1, 1e-15), clip = "on") +
  theme(strip.background = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA))

p_dndscv_leat <- tmp %>%
  filter(phenotype == "leat_α") %>%
  ggplot(aes(x = or_case_ctrl, y = qglobal_cv, color = significance)) +
  geom_point() +
  geom_hline(yintercept = label_threshold, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_label_repel(data = . %>% group_by(phenotype) %>% slice_min(qglobal_cv, n = 99, with_ties = F) %>% filter(or_case_ctrl > 0 & qglobal_cv < label_threshold), 
                   aes(label = sprintf('%s', gene_name), fontface = "italic"), 
                   hjust = 0,
                   label.size = .1,
                   min.segment.length = 0, 
                   label.padding = .1,
                   force = 1,
                   force_pull = .5,
                   size = 2,
                   show.legend = FALSE) +
  scale_y_continuous(trans = reverselog_trans(10), 
                     breaks = c(1e0, label_threshold, 1e-5, 1e-10, 1e-15),
                     oob = scales::oob_squish_infinite,
                     labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(trans = "log10", 
                     oob = scales::oob_squish_infinite,
                     labels = function(x) ifelse(x == 0, "0", x)) +
  scale_colour_identity() +
  theme_classic(base_size = 7) +
  xlab("OR Case vs. Control") +
  ylab("dNdScv Q-value") +
  coord_cartesian(xlim = c(-10, 100), ylim = c(1, 1e-15), clip = "on") +
  theme(strip.background = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA))

p_dndscv_subset <- cowplot::plot_grid(p_dndscv_mcd, 
                                      p_dndscv_leat, 
                                      nrow = 2,
                                      align = "hv")
