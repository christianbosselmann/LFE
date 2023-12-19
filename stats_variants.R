## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## Variant-level statistics of dNdScv meta-analysis calls versus population
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
                 readxl)

### DATA -----------------------------------------------------------------------
## ANNOVAR-annotated meta-analysis calls
df <- read_csv("data/MASTERLIST_DNDSCV_CALLS_SIMPLIFIED_2023-08-09.csv")

# remove IGV-confirmed seq artifact from calls
df <- df %>%
  filter(!cosmic_aa %like% "H304Y")

### ANALYSIS -------------------------------------------------------------------
# define gene of interest
gene <- "EGFR"

# split into sub-dataframes for gene of interest
df_subset <- df %>%
  filter(Gene.refGene == gene)

# fix col name
df_subset <- df_subset %>%
  rename(CADD_PHRED = CADD13_PHRED)

# get annotated variants from gnomAD and ClinVar
df_gnomad <- read_tsv(paste0("data/clinvar_gnomad/gnomAD_", gene, "_annotated.txt"))
df_clinvar_blb <- read_tsv(paste0("data/clinvar_gnomad/ClinVar_", gene, "_annotated.txt")) %>%
  filter(ClinicalSignificance %ilike% "benign")
df_clinvar_plp <- read_tsv(paste0("data/clinvar_gnomad/ClinVar_", gene, "_annotated.txt")) %>%
  filter(ClinicalSignificance %ilike% "pathogenic")

# merge
scores <- c("MTR", "EVE", "REVEL", "CADD_PHRED")
df_analysis <- rbind(tibble(group = "This study", df_subset[ ,scores]),
                     tibble(group = "gnomAD", df_gnomad[ ,scores]),
                     tibble(group = "ClinVar\nbenign", df_clinvar_blb[ ,scores]),
                     tibble(group = "ClinVar\npathogenic", df_clinvar_plp[ ,scores]))

# fix factor levels and pivot to long 
df_analysis <- df_analysis %>%
  mutate(group = factor(group, levels = c("This study", "ClinVar\npathogenic", 
                                          "ClinVar\nbenign", "gnomAD"))) %>%
  pivot_longer(cols = c(MTR, EVE, REVEL, CADD_PHRED), names_to = "score") 

p <- df_analysis %>%
  ggplot(aes(x = group, y = value, color = group, fill = group)) + 
  ggrain::geom_rain(alpha = .5,
                    point.args = list(size = .3),
                    point.args.pos = list(position = position_nudge(x = -0.15)),
                    boxplot.args = list(fill = "white", outlier.shape = NA, width = .15),
                    boxplot.args.pos = list(position = position_nudge(x = 0.0)),
                    violin.args.pos = list(side = "r", width = 0.7, position = position_nudge(x = 0.15))) +
  geom_pwc(method = "wilcox_test",
           ref.group = "This study",
           label = "p.signif",
           size = 0.2,
           label.size = 2.5) +
  scale_y_continuous(expand = expand_scale(mult = c(.01, .1))) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  xlab("") +
  ylab("") +
  coord_cartesian(expand = TRUE) +
  facet_wrap(score ~ .,
             scales = "free")

# export
pdf(paste0("output/variant_stats_", gene, ".pdf"), height = 4, width = 4)
p
dev.off()

