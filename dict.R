## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## Dictionary of constants and strings
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
librarian::shelf(readxl,
                 tidyverse,
                 janitor)

### DICT -----------------------------------------------------------------------
# established LFE associated genes
df_known <- readxl::read_excel("data/known_LFE_genes.xlsx") %>%
  janitor::clean_names() %>%
  select(gene_1, pheno) %>%
  rename(gene = gene_1)

vec_known <- df_known %>%
  pull(gene) %>%
  unique

# plus candidates from dNdScv
vec_known_and_hits <- c(vec_known, "DYRK1A", "EGFR")

# manually curated list of established genes with prior statistical support,
# established genes with first-time statistical support in this study, and
# novel genes
tmp <- readxl::read_excel("data/dndscv_significant_subset.xlsx", sheet = 2)
df_classification <- tmp %>% 
  distinct(Gene, Category)
