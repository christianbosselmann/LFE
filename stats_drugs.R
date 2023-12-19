## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## This script is used to support therapeutic actionability of the targets
## in this study by their number of ongoing trials and drug assocations
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
                 readxl)

### CLINICAL TRIALS ------------------------------------------------------------
# data: ClinicalTrials.gov full export for each target
paths <- list.files("data/ClinicalTrials", full.names = TRUE)
names <- basename(paths) %>% gsub(pattern = ".csv", "", .)
ls_ct <- lapply(paths, read_csv)
names(ls_ct) <- names
df_ct <- rbindlist(ls_ct, idcol = "target")

# sort for interventional trials not yet completed
df_ct <- df_ct %>%
  filter(`Study Type` == "INTERVENTIONAL") %>%
  filter(`Study Status` %like% "RECRUIT")

# summarize number of trials and each NCT number for each target
df_ct_res <- df_ct %>%
  group_by(target) %>%
  summarize(number_of_trials = n(),
            trial_ids = paste0(`NCT Number`, collapse = ";"))

# plot
p_ct <- df_ct_res %>%
  ggplot(aes(x = reorder(target, -number_of_trials), y = number_of_trials)) +
  geom_col() +
  geom_text(aes(label = number_of_trials, vjust = -.5), size = 2) +
  coord_cartesian() +
  scale_y_continuous(expand = expansion(add = c(0, 200))) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(face = "italic")) +
  guides(x = guide_axis(angle = 45)) +
  xlab("Target") +
  ylab("Ongoing trials on ClinicialTrials.gov")

### TARGET-DRUG ASSOCATIONS: DRUGCENTRAL ---------------------------------------
# get database
df_dc <- read_tsv("data/DrugCentral/drug.target.interaction.tsv")

# names of FDA approved medication
vec_fda <- read_csv("data/DrugCentral/FDA_Approved.csv", col_names = F) %>%
  pull(X2)

# subset DrugCentral database to our targets
df_dc <- df_dc %>%
  filter(GENE %in% names)

# annotate if drug is FDA approved for any indication
df_dc <- df_dc %>%
  mutate(fda_approved = ifelse(DRUG_NAME %in% vec_fda, 1, 0))

# get summary statistics: how many drugs per target, how are they called, are they FDA approved?
df_dc_res <- df_dc %>%
  group_by(GENE) %>%
  summarize(number_of_drugs = n(),
            number_of_fda_drugs = sum(fda_approved),
            drug_names = paste0(DRUG_NAME, collapse = ";"))

# plot
p_dc <- df_dc_res %>%
  mutate(non_fda_drugs = number_of_drugs - number_of_fda_drugs) %>%
  pivot_longer(cols = c(number_of_fda_drugs, non_fda_drugs)) %>%
  filter(value > 0) %>%
  ggplot(aes(x = reorder(GENE, -number_of_drugs), y = value, fill = name)) +
  geom_col(position = "stack") +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 2) +
  coord_cartesian() +
  scale_y_continuous(expand = expansion(add = c(0, 5))) +
  scale_fill_brewer("FDA approval status", 
                    palette = "Set2", 
                    direction = -1,
                    labels = c("Not approved", "Approved")) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(face = "italic")) +
  guides(x = guide_axis(angle = 45)) +
  xlab("Target") +
  ylab("Target-drug associations")

### EXPORT PLOTS ---------------------------------------------------------------
tmp <- cowplot::plot_grid(p_ct, p_dc, 
                          nrow = 1,
                          rel_widths = c(.45, .55))

pdf("output/drugs.pdf", height = 3, width = 6)
tmp
dev.off()
