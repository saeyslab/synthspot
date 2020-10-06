library(Seurat)
library(tidyverse)

seurat_obj= readRDS("../synthetic_data/seuratObj_subsampled.rds")
usethis::use_data(seurat_obj,overwrite = T, compress = "bzip2")




