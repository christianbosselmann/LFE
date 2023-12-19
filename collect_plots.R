## -----------------------------------------------------------------------------
##
## Somatic variant detection in lesional focal epilepsy
## This script collects plots and finalizes them for publication
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
                 thematic,
                 ggrepel,
                 scales,
                 cowplot)

# seed
set.seed(42)

## plot dimensions
nature_width_full <- 7.08661 # 180 mm
nature_height_short <- 7.28346 # 185 mm
nature_height_medium <- 8.26772 # 210 mm
nature_height_long <- 8.85827 # 225 mm

### FIGURE 1 -------------------------------------------------------------------
source("plot_dndscv.R")
source("plot_dndscv_qq.R")
source("plot_dndscv_scatter.R")

p_figure_1 <- cowplot::plot_grid(p_dndscv_subset, p_dndscv_qq, p_dndscv_scatter,
                   nrow = 1,
                   align = "hv")

pdf("output/Figure1.pdf", 
    width = nature_width_full-1, 
    height = nature_height_medium/2) 

p_figure_1

dev.off()
