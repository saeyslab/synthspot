<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

synthvisium
===========

<!-- badges: start -->

[![R build
status](https://github.com/browaeysrobin/synthvisium/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/browaeysrobin/synthvisium/actions)
[![Coverage
Status](https://codecov.io/gh/browaeysrobin/synthvisium/branch/master/graph/badge.svg?token=NKZBMJJDYA)](https://codecov.io/gh/browaeysrobin/synthvisium)
<!-- badges: end -->

**synthvisium: an R package to generate synthetic Visium-like spot data
based by generating artificial mixtures of cells from input scRNAseq
data.** The goal of this package is to give users the possibilities to
generate different types of synthetic data for evaluation of spot
deconvolution and label transfer methods.

Main functionalities of synthvisium
-----------------------------------

Add later.

Installation of synthvisium
---------------------------

Installation typically takes a few minutes, depending on the number of
dependencies that has already been installed on your pc.

You can install synthvisium (and required dependencies) from github
with:

    # install.packages("devtools")
    devtools::install_github("satijalab/seurat", ref = "spatial")
    devtools::install_github("saeyslab/nichenetr")
    devtools::install_github("browaeysrobin/synthvisium")

synthvisium was tested on both Windows and Linux (most recently tested R
version: R 4.0.2)

Learning to use synthvisium
---------------------------

In the following vignettes, you can find how to perform spatially aware
cell-cell communication inference based on integration of 10X Visium and
scRNAseq data:

-   [Basic computational analysis from spatial transcriptomics data:
    spatial regions approach](vignettes/basic_analysis.md):
    `vignette("basic_analysis", package="synthvisium")`

References
----------

Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular
communication by linking ligands to target genes. Nat Methods (2019)
<a href="doi:10.1038/s41592-019-0667-5" class="uri">doi:10.1038/s41592-019-0667-5</a>
