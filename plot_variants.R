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
                 readxl,
                 g3viz,
                 drawProteins)

# helper fns and dict
source("funs.R")
source("dict.R")

### DATA -----------------------------------------------------------------------
## ANNOVAR-annotated mega-analysis calls
df <- read_csv("data/MASTERLIST_DNDSCV_CALLS_SIMPLIFIED_2023-08-09.csv")

# split into sub-dataframes for genes of interest
df_dyrk1a <- df %>%
  filter(Gene.refGene %in% "DYRK1A")

df_egfr <- df %>%
  filter(Gene.refGene %in% "EGFR")

# extract MANE transcript amino acid exchange for each variant
# DYRK1A; ENST00000647188.2, contains ENSP00000494572 and matches to NM_001347721.2 and NP_001334650.1
df_dyrk1a <- df_dyrk1a %>%
  separate_wider_delim(AAChange.refGene, delim = ",", names_sep = "_", too_many = "debug", too_few = "debug") %>%
  rename(AAChange.refGene = AAChange.refGene_4) %>%
  select(-starts_with("AAChange.refGene_")) %>%
  mutate(AAChange = sub(".*p.", "", AAChange.refGene)) %>%
  mutate(pos = parse_number(AAChange))

# EGFR; ENST00000275493.7, contains ENSP00000275493 and matches to NM_005228.5 and NP_005219.2
df_egfr <- df_egfr %>%
  separate_wider_delim(AAChange.refGene, delim = ",", names_sep = "_", too_many = "debug", too_few = "debug") %>%
  rename(AAChange.refGene = AAChange.refGene_6) %>%
  select(-starts_with("AAChange.refGene_")) %>%
  mutate(AAChange = sub(".*p.", "", AAChange.refGene)) %>%
  mutate(pos = parse_number(AAChange)) %>%
  filter(!AAChange == "H304Y") # post-hoc exclusion of a sequencing artifact confirmed on IGV

## pfam domains of proteins of interest
# DYRK1A UniProt: Q13627
pfam_dyrk1a <- uniprot2pfam("Q13627")
# EGFR UniProt: P00533
pfam_egfr <- uniprot2pfam("P00533")

## ClinVar variant distributions
clinvar_dyrk1a <- read_csv("data/clinvar_DYRK1A.csv")
clinvar_dyrk1a <- clinvar_dyrk1a %>%
  filter(Type == "SNV") %>%
  mutate(pos = parse_number(HGVSp, locale = locale(grouping_mark = ". ", decimal_mark = ",")))

clinvar_egfr <- read_csv("data/clinvar_EGFR.csv")
clinvar_egfr <- clinvar_egfr %>%
  filter(Type == "SNV") %>%
  mutate(pos = parse_number(HGVSp, locale = locale(grouping_mark = ". ", decimal_mark = ",")))

## in silico variant distribution
insilico_dyrk1a <- read_csv("data/insilico_DYRK1A.csv") %>%
  mutate(pos = parse_number(HGVSp, locale = locale(grouping_mark = ". ", decimal_mark = ",")))

insilico_egfr <- read_csv("data/insilico_EGFR.csv") %>%
  mutate(pos = parse_number(HGVSp, locale = locale(grouping_mark = ". ", decimal_mark = ",")))

## MTR distribution
mtr_dyrk1a <- read_tsv("data/mtr_DYRK1a.csv")
mtr_egfr <- read_tsv("data/mtr_EGFR.csv")

## IUPRED distribution
iupred_dyrk1a <- read_table("data/iupred_DYRK1A.txt")[,1:3]
colnames(iupred_dyrk1a) <- c("pos", "aa", "iupred")

iupred_egfr <- read_table("data/iupred_EGFR.txt")[,1:3]
colnames(iupred_egfr) <- c("pos", "aa", "iupred")

### PLOT: DYRK1A ---------------------------------------------------------------
protein_length <- 755

# ClinVar density
p1 <- clinvar_dyrk1a %>%
  filter(!Significance %like% "Uncertain") %>%
  filter(!Significance %like% "Conflicting") %>%
  mutate(group = ifelse(Significance %ilike% "benign", "BLB", "PLP")) %>%
  distinct(pos, group) %>%
  na.omit %>%
  ggplot(aes(x = pos, group = group, fill = group, color = group)) +
  geom_density(alpha = 0.5, adjust = 1) +
  geom_rug(outside = TRUE, sides = "b", length = unit(0.2, "npc")) +
  coord_cartesian(xlim = c(0, protein_length), clip = "off", expand = FALSE) +
  theme_classic() +
  theme(legend.position = c(.9,.75),
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("ClinVar") +
  scale_fill_manual(values = c("lightgray", "#7d0013")) +
  scale_color_manual(values = c("lightgray", "#7d0013")) 
p1 <- makeSmallLegend(p1)

# gnomAD density
p2 <- insilico_dyrk1a %>%
  distinct(gnomAD_exomes_AC, pos) %>%
  na.omit %>%
  ggplot(aes(x = pos)) +
  geom_density(alpha = 0.5) +
  geom_rug(outside = TRUE, sides = "b", length = unit(0.2, "npc")) +
  coord_cartesian(xlim = c(0, protein_length), clip = "off", expand = FALSE) +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("gnomAD")

# MTR map
p3 <- mtr_dyrk1a %>%
  ggplot(aes(x = protein_position, y = mtr)) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = quantile(mtr_dyrk1a$mtr, .5), linetype = "dashed", color = "#7d0013") +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0.5, quantile(mtr_dyrk1a$mtr, .5)[[1]]),
                     labels = scales::label_number(accuracy = 0.1)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("MTR")

# IUPRED map
p4 <- iupred_dyrk1a %>%
  ggplot(aes(x = pos, y = iupred)) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#7d0013") +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0, 0.5), limits = c(0, max(iupred_dyrk1a$iupred))) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("IUPRED")

# lollipop plot
p5 <- df_dyrk1a %>%
  group_by(pos) %>%
  summarize(n = n()) %>%
  left_join(df_dyrk1a[,c("pos", "AAChange")]) %>%
  distinct %>%
  ggplot(aes(x = pos, y = n, label = AAChange)) +
  geom_point() + 
  geom_segment(aes(x = pos, xend = pos, y = 0, yend = n)) +
  ggrepel::geom_label_repel(size = 3, box.padding = 0.1, label.padding = .1) +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0, 4), limits = c(0, 4.5)) +
  theme_classic() +
  xlab("") +
  ylab("Variants") 

# protein plot
dp_dyrk1a <- drawProteins::get_features("Q13627")
dp_dyrk1a <- drawProteins::feature_to_dataframe(dp_dyrk1a)
pfam_dyrk1a <- draw_canvas(dp_dyrk1a)
pfam_dyrk1a <- draw_chains(pfam_dyrk1a, dp_dyrk1a, 
                           fill = "lightgray",
                           outline = "lightgray",
                           label_chains = FALSE)
pfam_dyrk1a <- draw_domains(pfam_dyrk1a, dp_dyrk1a,
                            label_domains = FALSE)
pfam_dyrk1a <- draw_phospho(pfam_dyrk1a, dp_dyrk1a)
pfam_dyrk1a <- pfam_dyrk1a + 
  theme_bw(base_size = 20) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.y.left = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(0, protein_length), ylim = c(0.65, 1.35), expand = FALSE) 
pfam_dyrk1a <- makeSmallLegend(pfam_dyrk1a)
pfam_dyrk1a_legend <- ggpubr::get_legend(pfam_dyrk1a)
pfam_dyrk1a_legend <- ggplotify::as.ggplot(pfam_dyrk1a_legend)
pfam_dyrk1a <- pfam_dyrk1a +
  theme(legend.position = "none")
pfam_dyrk1a

# bar charts for candidate variant scores
p6 <- df_dyrk1a %>%
  distinct(pos, REVEL) %>%
  ggplot(aes(x = pos, y = REVEL, fill = "#7d0013")) +
  geom_col(width = 7.5) +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_fill_manual(values = "#7d0013") +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("REVEL")

p7 <- df_dyrk1a %>%
  distinct(pos, EVE) %>%
  ggplot(aes(x = pos, y = EVE, fill = "#7d0013")) +
  geom_col(width = 7.5) +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_fill_manual(values = "#7d0013") +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("EVE")

p8 <- df_dyrk1a %>%
  distinct(pos, CADD13_PHRED) %>%
  ggplot(aes(x = pos, y = CADD13_PHRED, fill = "#7d0013")) +
  geom_col(width = 7.5) +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_fill_manual(values = "#7d0013") +
  scale_y_continuous(breaks = c(0, 20), limits = c(0, 35)) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("CADD")

# compose plot
plot_dyrk1a <- cowplot::plot_grid(p1 + theme(plot.margin = unit(c(0.1, 0.3, 0, 0), "cm")), 
                                  p2 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                  p3 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                  p4 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                  p6 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                  p7 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                  p8 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                  p5 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                  pfam_dyrk1a + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                  ncol = 1, 
                                  align = "hv")

pdf("output/plot_dyrk1a.pdf", height = 8*1.5, width = 4*1.5)
plot_dyrk1a
dev.off()

pdf("output/plot_dyrk1a_legend.pdf", height = 1*1.5, width = 3*1.5)
pfam_dyrk1a_legend
dev.off()

### PLOT: EGFR ---------------------------------------------------------------
protein_length <- 1211

# ClinVar density
p1 <- clinvar_egfr %>%
  filter(!Significance %like% "Uncertain") %>%
  filter(!Significance %like% "Conflicting") %>%
  mutate(group = ifelse(Significance %ilike% "benign", "BLB", "PLP")) %>%
  distinct(pos, group) %>%
  na.omit %>%
  ggplot(aes(x = pos, group = group, fill = group, color = group)) +
  geom_density(alpha = 0.5, adjust = 1) +
  geom_rug(outside = TRUE, sides = "b", length = unit(0.2, "npc")) +
  coord_cartesian(xlim = c(0, protein_length), clip = "off", expand = FALSE) +
  theme_classic() +
  theme(legend.position = c(.9,.75),
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("ClinVar") +
  scale_fill_manual(values = c("lightgray", "#7d0013")) +
  scale_color_manual(values = c("lightgray", "#7d0013")) 
p1 <- makeSmallLegend(p1)

# gnomAD density
p2 <- insilico_egfr %>%
  distinct(gnomAD_exomes_AC, pos) %>%
  na.omit %>%
  ggplot(aes(x = pos)) +
  geom_density(alpha = 0.5) +
  geom_rug(outside = TRUE, sides = "b", length = unit(0.2, "npc")) +
  coord_cartesian(xlim = c(0, protein_length), clip = "off", expand = FALSE) +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("gnomAD")

# MTR map
p3 <- mtr_egfr %>%
  ggplot(aes(x = protein_position, y = mtr)) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = quantile(mtr_egfr$mtr, .5), linetype = "dashed", color = "#7d0013") +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0.5, quantile(mtr_egfr$mtr, .5)[[1]]),
                     labels = scales::label_number(accuracy = 0.1)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("MTR")

# IUPRED map
p4 <- iupred_egfr %>%
  ggplot(aes(x = pos, y = iupred)) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#7d0013") +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0, 0.5), limits = c(0, max(iupred_egfr$iupred))) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("IUPRED")

# lollipop plot
p5 <- df_egfr  %>%
  group_by(pos) %>%
  summarize(n = n()) %>%
  left_join(df_egfr[,c("pos", "AAChange")]) %>%
  distinct %>%
  ggplot(aes(x = pos, y = n, label = AAChange)) +
  geom_point() + 
  geom_segment(aes(x = pos, xend = pos, y = 0, yend = n)) +
  ggrepel::geom_label_repel(size = 3, box.padding = 0.1, label.padding = .1) +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0, 4), limits = c(0, 4.5)) +
  theme_classic() +
  xlab("") +
  ylab("Variants") 

# protein plot
dp_egfr <- drawProteins::get_features("P00533")
dp_egfr <- drawProteins::feature_to_dataframe(dp_egfr)
pfam_egfr <- draw_canvas(dp_egfr)
pfam_egfr <- draw_chains(pfam_egfr, dp_egfr, 
                         fill = "lightgray",
                         outline = "lightgray",
                         label_chains = FALSE)
pfam_egfr <- draw_domains(pfam_egfr, dp_egfr,
                          label_domains = FALSE)
pfam_egfr <- draw_repeat(pfam_egfr, dp_egfr,
                         label_repeats = FALSE)
pfam_egfr <- draw_regions(pfam_egfr, dp_egfr)
pfam_egfr <- draw_phospho(pfam_egfr, dp_egfr)
pfam_egfr <- pfam_egfr + 
  theme_bw(base_size = 20) +
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(),
        axis.title.y.left = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(0, protein_length), ylim = c(0.65, 1.35), expand = FALSE) 
pfam_egfr <- makeSmallLegend(pfam_egfr)
pfam_egfr_legend <- ggpubr::get_legend(pfam_egfr)
pfam_egfr_legend <- ggplotify::as.ggplot(pfam_egfr_legend)
pfam_egfr <- pfam_egfr +
  theme(legend.position = "none")

# bar charts for candidate variant scores
p6 <- df_egfr %>%
  distinct(pos, REVEL) %>%
  ggplot(aes(x = pos, y = REVEL, fill = "#7d0013")) +
  geom_col(width = 7.5) +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
  scale_fill_manual(values = "#7d0013") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("REVEL")

p7 <- df_egfr %>%
  distinct(pos, EVE) %>%
  ggplot(aes(x = pos, y = EVE, fill = "#7d0013")) +
  geom_col(width = 7.5) +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
  scale_fill_manual(values = "#7d0013") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("EVE")

p8 <- df_egfr %>%
  distinct(pos, CADD13_PHRED) %>%
  ggplot(aes(x = pos, y = CADD13_PHRED, fill = "#7d0013")) +
  geom_col(width = 7.5) +
  coord_cartesian(xlim = c(0, protein_length), expand = FALSE) +
  scale_y_continuous(breaks = c(0, 20), limits = c(0, 35)) +
  scale_fill_manual(values = "#7d0013") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("CADD")

# compose plot
plot_egfr <- cowplot::plot_grid(p1 + theme(plot.margin = unit(c(0.1, 0.3, 0, 0), "cm")), # trbl 
                                p2 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                p3 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                p4 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                p6 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                p7 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                p8 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                p5 + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                pfam_egfr + theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm")), 
                                ncol = 1, 
                                align = "hv")

pdf("output/plot_egfr.pdf", height = 8*1.5, width = 4*1.5)
plot_egfr
dev.off()

pdf("output/plot_egfr_legend.pdf", height = 1*1.5, width = 3*1.5)
pfam_egfr_legend
dev.off()

### EXPORT DATA ----------------------------------------------------------------
## combine candidate dataframes for structural analysis
df_export <- rbind(df_dyrk1a, df_egfr)

write_csv(df_export, "data/variants_in_candidate_genes.csv")
