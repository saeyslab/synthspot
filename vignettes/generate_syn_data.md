Generating Visium-like spot data
================
Robin Browaeys
2020-10-08

<!-- github markdown built using 
rmarkdown::render("vignettes/generate_syn_data.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to use synthvisium to generate
synthetic Visium-like spot data. First we will demonstrate it by using a
toy scRNAseq dataset as input (internal Seurat object of downsampled
brain cortex data). Later on, we will demonstrate it on a real scRNAseq
dataset (brain cortex) that can be downloaded from Zenodo. Using this
dataset as input for your own benchmark is recommended.

``` r
library(Seurat)
library(synthvisium)
library(dplyr)
```

## Generating Visium-like spot data from toy input scRNAseq data.

#### Explore input scRNAseq data

First we will explore the toy scRNAseq dataset (`seurat_obj`) and show
where the cell type information is stored in this object (`subclass`
column of the meta data)

``` r
# scrnaseq data
DimPlot(seurat_obj, group.by = "subclass", label = T)
```

![](generate_syn_data_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

#### Generate the synthetic data

We provide one function to generate a dataset and its meta data. As
example here, we will generate a dataset consisting of 5 different
regions with each between 5 or 20 spots (and on average 30000 counts in
total per spot). The dataset type that will be created here is called
`artificial_diverse_distinct`, and is one of the 13 different
implemented options. As can be read in the documentation of the function
(`?generate_synthetic_visium`), here each artifical prior region has a
distinct set of cell types, and cell type frequencies are randomly
different for each cell type in each region.

``` r
synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_diverse_distinct", clust_var = "subclass", n_regions = 5, n_spots_min = 5, n_spots_max = 20, visium_mean = 30000, visium_sd = 8000)
```

#### Explore the synthetic data object

Now we will go over the different components of this object:

First, we have the gene expression counts for each spot in the `counts`
sublist (genes in rows, spots in columns):

``` r
synthetic_visium_data$counts %>% as.matrix() %>% .[1:5,1:5]
##      priorregion1_spot_1 priorregion1_spot_2 priorregion1_spot_3 priorregion1_spot_4 priorregion1_spot_5
## Vip                    0                   4                   2                   1                4242
## Sst                    1                   0                   0                   1                   1
## Npy                 1908                3339                9498                3696                2574
## Tac2                   0                 416                   0                  83                1063
## Crh                    0                   0                 588                   0                 170
```

Then, you can also see the cell type composition (absolute numbers of
cells and relative proportions) of each generated spot in the
`spot_composition` `relative_spot_composition` sublists (names and prior
regions of the spot in the `name` and `region` columns):

``` r
synthetic_visium_data$spot_composition %>% .[1:10,]
##    L2/3.IT L4 L5.IT L5.PT L6.CT L6.IT L6b Lamp5 NP Pvalb Sst Vip                name       region
## 1        0  0     0     0     0     0   0     7  0     0   0   0 priorregion1_spot_1 priorregion1
## 2        0  0     0     0     0     0   0     6  0     0   0   0 priorregion1_spot_2 priorregion1
## 3        0  0     0     0     0     0   0     4  0     0   0   1 priorregion1_spot_3 priorregion1
## 4        0  0     0     0     0     0   0     6  0     0   0   0 priorregion1_spot_4 priorregion1
## 5        0  0     0     0     0     0   0     4  0     0   0   1 priorregion1_spot_5 priorregion1
## 6        0  0     0     0     0     0   0     4  0     0   0   1 priorregion1_spot_6 priorregion1
## 7        0  0     0     0     0     0   0     7  0     0   0   0 priorregion1_spot_7 priorregion1
## 8        0  0     2     0     0     0   0     0  0     3   0   0 priorregion2_spot_1 priorregion2
## 9        0  0     0     0     0     0   0     0  0     3   0   0 priorregion2_spot_2 priorregion2
## 10       0  0     3     0     0     0   0     0  0     5   0   0 priorregion2_spot_3 priorregion2
synthetic_visium_data$relative_spot_composition %>% .[1:10,]
##    L2/3.IT L4 L5.IT L5.PT L6.CT L6.IT L6b Lamp5 NP Pvalb Sst Vip                name       region
## 1        0  0 0.000     0     0     0   0   1.0  0 0.000   0 0.0 priorregion1_spot_1 priorregion1
## 2        0  0 0.000     0     0     0   0   1.0  0 0.000   0 0.0 priorregion1_spot_2 priorregion1
## 3        0  0 0.000     0     0     0   0   0.8  0 0.000   0 0.2 priorregion1_spot_3 priorregion1
## 4        0  0 0.000     0     0     0   0   1.0  0 0.000   0 0.0 priorregion1_spot_4 priorregion1
## 5        0  0 0.000     0     0     0   0   0.8  0 0.000   0 0.2 priorregion1_spot_5 priorregion1
## 6        0  0 0.000     0     0     0   0   0.8  0 0.000   0 0.2 priorregion1_spot_6 priorregion1
## 7        0  0 0.000     0     0     0   0   1.0  0 0.000   0 0.0 priorregion1_spot_7 priorregion1
## 8        0  0 0.400     0     0     0   0   0.0  0 0.600   0 0.0 priorregion2_spot_1 priorregion2
## 9        0  0 0.000     0     0     0   0   0.0  0 1.000   0 0.0 priorregion2_spot_2 priorregion2
## 10       0  0 0.375     0     0     0   0   0.0  0 0.625   0 0.0 priorregion2_spot_3 priorregion2
```

If you want to see which cell type is present in which region, and with
which prior frequency, you can look at the `gold_standard_priorregion`
sublist.

``` r
synthetic_visium_data$gold_standard_priorregion %>% head()
## # A tibble: 6 x 4
##   prior_region celltype  freq present
##   <chr>        <chr>    <dbl> <lgl>  
## 1 priorregion1 Vip        0.1 TRUE   
## 2 priorregion1 Lamp5      0.9 TRUE   
## 3 priorregion1 Sst        0   FALSE  
## 4 priorregion1 L4         0   FALSE  
## 5 priorregion1 L5 IT      0   FALSE  
## 6 priorregion1 L6 IT      0   FALSE
```

We also provide a data frame with properties of this specific dataset.
This might be useful in a benchmark to check the influence of these
properties on performance.

``` r
synthetic_visium_data$dataset_properties 
## # A tibble: 1 x 8
##   dataset_id                   dataset_type                real_artificial uniform_diverse distinct_overlap dominant_celltype missing_celltype real_region_var
##   <chr>                        <chr>                       <chr>           <chr>           <chr>            <lgl>             <lgl>            <lgl>          
## 1 artificial_diverse_distinct1 artificial_diverse_distinct artificial      diverse         distinct         FALSE             FALSE            NA
```

## Generating Visium-like spot data from toy input scRNAseq data that contains information of the anatomical regions of cells.

Just pooling cells from different cell types as done in the previous
example could have the (theoretical) disadvantage that spatial
information inherently encoded in the scRNAseq data is lost. This is
because cells extracted from a similar anatomical region could be
influenced by the same microenvironmental signal, resulting in a shared
spatial gene expression pattern for cells in the same region, which
could be reflected in the input scRNAseq data. In case you know the
anatomical origin from cells in the input scRNAseq data, it could be
useful to generate synthetic Visium-spot data by only sampling cells in
one spot if those cells belong to the same anatomical region.

In the brain cortex dataset, we have this information stored in the meta
data (`brain_subregion` column)

``` r
# scrnaseq data
DimPlot(seurat_obj, group.by = "brain_subregion")
```

![](generate_syn_data_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

As the second example, we will now generate a dataset consisting of
artificial “L1 - L2/3 - L4 - L5 - L6” regions with each between 5 or 20
spots (and on average 30000 counts in total per spot). The dataset type
that will be created here is called `real`, and is one of the 13
different implemented options. As can be read in the documentation of
the function (`?generate_synthetic_visium`), this dataset type will
exploit the anatomical information and cell type frequencies per region
will be based on the number of cells in the corresponding region in the
input scRNAseq data.

``` r
synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real", clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
```

## Generating Visium-like spot data from real input scRNAseq data.

The real input scRNAseq (Allen Brain Atlas brain cortex) can be
downloaded from Zenodo (<https://zenodo.org/record/4072975>).

``` r
# seurat_obj = readRDS(url("https://zenodo.org/record/4072975/files/seurat_obj_scrnaseq_cortex_filtered.rds")) # direct read from URL - will take a long time
# seurat_obj = readRDS(path/seurat_obj_scrnaseq_cortex_filtered.rds")) # after downloading and saving locally
```

Now, we will create a synthetic dataset based on this real dataset and
with a more realistic number of spots. This will take some time to run.
As another example dataset type, we pick here `real_top2_overlap`.

``` r
synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real_top2_overlap", clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 50, n_spots_max = 200, visium_mean = 20000, visium_sd = 5000)
```

Check the following vignette to find a short demonstration of a basic
analysis you can do on this synthetic Visium data. We also show how to
do the label transfer and compare gold standard labels with predicted
labels. [Analyzing synthetic Visium-like spot data for evaluating label
transfer methods](analyze_syn_data.md):`vignette("analyze_syn_data",
package="synthvisium")`.

## Information about the different dataset types.

Currently, there are 13 possible dataset types implemented:

``` r
possible_dataset_types = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform", "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",  "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium")
```

Explanation of each of the types:

  - **real**: for each region (defined in scRNAseq data - region\_var -
    and thus **real** ), cell type frequencies are based on the number
    of cells in the corresponding region in the input scRNAseq data
  - **real\_top1**: for each region (defined in scRNAseq data and thus
    **real** ), we will only consider cell types that are most
    frequently present in that 1 region ( **top 1** ). Cell type
    frequencies are based on the number of cells in the corresponding
    region in the input scRNAseq data
  - **real\_top1\_uniform**: for each region (defined in scRNAseq data -
    region\_var - and thus **real** ), we will only consider cell types
    that are most frequently present in that 1 region ( **top 1** ).
    Cell type frequencies are equal for each cell type assigned to a
    region ( **uniform** )
  - **real\_top2\_overlap**: for each region (defined in scRNAseq data -
    region\_var - and thus **real** ), we will only consider cell types
    that if this particular region is among the 2 regions in which they
    are most abundant ( **top 2** ). As a consequence, there will be
    some **overlap** in cell type presence between two regions. Cell
    type frequencies are based on the number of cells in the
    corresponding region in the input scRNAseq data
  - **real\_top2\_overlap\_uniform**: for each region (defined in
    scRNAseq data - region\_var - and thus **real** ), we will only
    consider cell types that if this particular region is among the 2
    regions in which they are most abundant ( **top 2** ). As a
    consequence, there will be some **overlap** in cell type presence
    between two regions. Cell type frequencies are equal for each cell
    type assigned to a region ( **uniform** )
  - **real\_missing\_celltypes\_visium**: the same as the **real**
    dataset type, but now some random input scRNAseq cell types will be
    left out, and thus not used to generate synthetic data.
  - **artificial\_uniform\_distinct**: celltype-prior region assignments
    are entirely **artificial** and thus not based on possibly present
    regional information in the input scRNAseq data. Each artifical
    prior region has a distinct set of cell types here, and cell type
    frequencies are the same for each cell type in each region.
  - **artificial\_diverse\_distinct**: celltype-prior region assignments
    are entirely **artificial** and thus not based on possibly present
    regional information in the input scRNAseq data. Each artifical
    prior region has a distinct set of cell types here, and cell type
    frequencies are different for each cell type in each region.  
  - **artificial\_uniform\_overlap**: celltype-prior region assignments
    are entirely **artificial** and thus not based on possibly present
    regional information in the input scRNAseq data. Each artifical
    prior region could have some overlap in cell type contribution with
    other regions, and cell type frequencies are the same for each cell
    type in each region.  
  - **artificial\_diverse\_overlap**: celltype-prior region assignments
    are entirely **artificial** and thus not based on possibly present
    regional information in the input scRNAseq data. Each artifical
    prior region could have some overlap in cell type contribution with
    other regions, and cell type frequencies are different for each cell
    type in each region.  
  - **artificial\_dominant\_celltype\_diverse**: celltype-prior region
    assignments are entirely **artificial** and thus not based on
    possibly present regional information in the input scRNAseq data.
    Each artifical prior region could have some overlap in cell type
    contribution with other regions, and cell type frequencies are
    different for each cell type in each region. There is one dominant
    celltype that is present in each region, and is much more abundant
    than other cell types (5-15x more).
  - **artificial\_partially\_dominant\_celltype\_diverse**:
    celltype-prior region assignments are entirely **artificial** and
    thus not based on possibly present regional information in the input
    scRNAseq data. Each artifical prior region could have some overlap
    in cell type contribution with other regions, and cell type
    frequencies are different for each cell type in each region. There
    is one dominant celltype that is present in all but one regions
    (partially dominant), and is much more abundant than other cell
    types, except for one region where it is equally abundant.
  - **artificial\_missing\_celltypes\_visium**: the same as the
    ‘artificial\_diverse\_overlap’ dataset type, but now some random
    input scRNAseq cell types will be left out, and thus not used to
    generate synthetic data.
