## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## This script plots the meta-analysis cohort size and phenotypes
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
df_calls <- read_csv("data/MASTERLIST_DNDSCV_CALLS_SIMPLIFIED_2023-08-09.csv")

### BAR PLOT -------------------------------------------------------------------
# cohort size and number of calls by data source for dNdScv
df_bar <- df_calls %>%
  group_by(Source) %>%
  summarise(n_sample = n_distinct(ID), n_calls = n())

p_bar <- df_bar %>%
  gather(n_sample, n_calls, -Source) %>%
  ggplot(aes(x = Source, y = n_calls, fill = n_sample, color = n_sample)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(trans = "log10", breaks = c(1, 100, 1000))+
  scale_x_discrete(labels = c("Baulac et al.\n(WES)", "Chung et al.\n(WES)",
                              "Panel cohort\n", "WES cohort\n"),
                   guide = guide_axis(angle = 30)) +
  scale_fill_brewer(labels = c("calls", "samples"), "Number of",
                    palette = "Greys") +
  scale_color_manual(values = c("black", "black")) +
  theme_classic(base_size = 7) +
  theme(legend.position = "right") +
  guides(color = "none") +
  coord_cartesian(expand = FALSE) +
  ylab("n") +
  xlab("")

pdf("output/bar.pdf", height = 2, width = 3)
p_bar
dev.off()

### PIE: PHENOTYPES ------------------------------------------------------------
# pie chart distribution of phenotype categories by source
df_pie_samples <- df_calls %>%
  group_by(Source, Main_pathology) %>%
  summarise(n = n_distinct(ID)) %>%
  mutate(Main_pathology = ifelse(Main_pathology %in% c("LEAT", "MCD"), Main_pathology, "Other")) %>%
  group_by(Source, Main_pathology) %>%
  mutate(n = sum(n)) %>%
  distinct

# fill up to ratio
df_pie_samples <- df_pie_samples %>%
  group_by(Source) %>%
  mutate(ratio = n/sum(n))

pie_names <- list(
  'Baulac_WES' = "Baulac et al.",
  'Chung_WES' = "Chung et al.",
  'Lal_panel' = "Panel cohort",
  'Lal_WES' = "WES cohort"
)

pie_labeller <- function(variable,value){
  return(pie_names[value])
}

p_pie <- df_pie_samples %>%
  ggplot(aes(x="", y = ratio, fill = Main_pathology)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Greys", "Pathology") +
  theme_void(base_size = 7) +
  facet_wrap(Source ~ ., ncol = 2,
             labeller = pie_labeller,
             strip.position = "bottom") +
  theme(strip.text.x.bottom =  element_text(angle = 0, vjust = 1.1, size = 7))

pdf("output/pie.pdf", height = 3, width = 3)
p_pie
dev.off()

### PIE: GENES -----------------------------------------------------------------
# pie chart distribution of phenotype categories by source
df_pie_genes <- df_calls %>%
  mutate(is_known = ifelse(Gene.refGene %in% vec_known, "Established genes", 
                           ifelse(Gene.refGene %in% c("DYRK1A", "EGFR"), "Candidate genes", 
                                  "Not associated"))) %>%
  group_by(is_known) %>%
  summarize(n = n())

p_pie_genes <- df_pie_genes %>%
  mutate(is_known = factor(is_known, levels = c("Not associated", "Established genes", "Candidate genes"))) %>%
  ggplot(aes(x = "", y = n, fill = is_known)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Greys", "Gene-disease association") +
  theme_void(base_size = 7) 

pdf("output/pie_genes.pdf", height = 3, width = 3)
p_pie_genes
dev.off()

