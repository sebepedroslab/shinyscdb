This package contains functions for single-cell atlas database shiny apps, including single species metacell atlases and comparative anlyisis between pairs of species.

To use these functions, each app needs to have `config.yaml` in which input files are specified, and `global.R` (example below), with libraries and global variables, including list of species (the species acronyms should match those used in config file).

```
# packages
options(repos = BiocManager::repositories())
library(BiocManager)
library(shiny)
library(data.table)
library(dplyr)
library(reshape2)
library(stringr)
library(glue)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(gplots)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(zoo)
library(kableExtra)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(patchwork)
library(rtracklayer)
library(shinyscdb)
library(shinycssloaders)

# config file
conf <- yaml::yaml.load_file("config.yaml", eval.expr=TRUE)

# species dictionary with '"full name" = acronym' pairs
sps <- c(
  "Hoilungia hongkongensis" = "Hhon",
  "Trichoplax adhaerens" = "Tadh",
  "Nematostella vectensis" = "Nvec"
)

# graphical params
heatmap_colors <- c("white","gray99","orange","orangered2","#520c52")

# env interactive heatmap
shiny_env = new.env()
```

