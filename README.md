<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

# synthspot

<!-- badges: start -->

[![R build
status](https://github.com/browaeysrobin/synthvisium/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/browaeysrobin/synthvisium/actions)
[![Coverage
Status](https://codecov.io/gh/browaeysrobin/synthvisium/branch/master/graph/badge.svg?token=MIOIRFJDQN)](https://codecov.io/gh/browaeysrobin/synthvisium)
<!-- badges: end -->

**synthspot: an R package to generate synthetic Visium-like spot data
by generating artificial mixtures of cells from input scRNAseq data.**
The goal of this package is to give users the possibilities to generate
different types of synthetic data for evaluation of spot deconvolution
and label transfer methods. Since these synthetic datasets just consist
of spots/mixtures of cells without modeling spatial aspects, they cannot
be used to evaluate methods that require a spatial component in the
data.

## How does this synthetic data generator work?

The goal of this simulator is to create a dataset that provide gene
expression information for ‘mini-bulk’ spots, similar to data generated
by the classic Spatial Transcriptomics and 10X Visium technology.
Because current spatial transcriptomics profile gene expression of spots
that typically contain between 2 and 10 cells, this simulator will
generate synthetic spot data by considering the gene expression of one
spot as a mixture of the expression of n cells (n between 2 and 10).

So to generate one spot, the simulator will sample n cells from input
scRNAseq data, sum their counts gene by gene, and downsample to a
‘standard Visium’ number of counts (because by summing counts of
multiple cells without downsampling, we would otherwise get
unrealistically high number of counts). n is currently a randomly picked
number between 2 an 10 (randomly sampled based on uniform distribution).
The target number of counts for downsampling is currently picked from a
normal distribution with mean and standard deviation based on real
Visium datasets (these values can be changed by the user).

Importantly, the sampling of cells from input scRNAseq data depends on
***a priori*** defined cell type frequencies. This means that cells from
scRNAseq-derived cell types with higher defined frequencies will be more
likely to be picked than cells from less frequent cell types. To make
more realistic datasets, the simulator allows for generating different
‘regions’, which are groups of spots with similar cell type frequencies.
By assigning different cell type frequencies to different ‘regions’,
label transfer methods to predict Visium-based region annotations of
single-cells can be applied and evaluated. In this simulator, several
‘dataset types’ with each different ways to define the prior cell type
frequencies are build in to represent different possible tissue types.
Details about these different types are given in the vignette.

Current limitations of this generator are that potential technical
differences between scRNAseq and Visium data are not yet incorporated
and that a region is just a group of spots and that a spatial
information layer is lacking.

## Installation of synthspot

Installation depends on the number of dependencies that has already been
installed on your pc, but should not take much time.

You can install synthspot (and required dependencies) from github
with:

    # install.packages("devtools")
    # install.packages("igraph", type = "binary") ### on Windows only, if problem with installing igraph
    devtools::install_github("saeyslab/synthspot")

synthspot has been tested on Windows, Mac and Linux (most recently
tested R version: R 4.0.2)

## Learning to use synthspot

In the following vignette, you can learn how to generate the synthetic
data by synthspot:

-   [Generating Visium-like spot data](vignettes/generate_syn_data.md):
    `vignette("generate_syn_data", package="synthspot")`

In the next vignette, you can find a short demonstration of a basic
analysis you can do on this synthetic Visium data. We also show how to
do the label transfer and compare gold standard labels with predicted
labels:

-   [Analyzing synthetic Visium-like spot data for evaluating label
    transfer methods](vignettes/analyze_syn_data.md):
    `vignette("analyze_syn_data", package="synthspot")`
