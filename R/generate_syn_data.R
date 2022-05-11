#' @title Prepare the generation of synthetic Visium/ST data by using a scRNAseq Seurat object as input and defining a priori defined cell type frequencies per region (group of spots).
#'
#' @description \code{make_region_celltype_assignment} will start from input scRNAseq data (Seurat object) and generate cell type frequencies that are used as base for generating synthetic Visium.
#' Via this functions you can generate different 'tissue types' by changing the way the cell type frequencies are defined. Read the argument description of `dataset_type` for seeing the different tissue type options.
#' If your scRNAseq data object already provides information about the region/location a cell was sampled from, you could use this information as well to generate synthetic Visium data that will only pool cells together if they originate from the same region.
#' @usage
#' make_region_celltype_assignment(seurat_obj, clust_var, n_regions, dataset_type, region_var = NULL, dataset_id = "1", n_spots_min = 50, n_spots_max = 500)
#'
#' @param seurat_obj The input scRNAseq data stored as Seurat object.
#' @param clust_var Name of the meta data column in which the cell type label of interest is defined.
#' @param n_regions Number of regions you want to define. If applicable: this can be NULL if you want to use the same number of regions as defined in a meta data column indicating the region of interest.
#' @param region_var Default NULL. If your scRNAseq data object already provides information about the region/location a cell was sampled from, you could change this argument to the meta data column indicating the region of interest.
#' @param dataset_id Character vector to define how you want to identify the dataset. The final identifier/name of the dataset will be the description of the dataset type followed by this character vector.
#' @param n_spots_min The minimum number of spots allowed per region. Default 50.
#' @param n_spots_max The maximum number of spots allowed per region. Default: 500.
#' @param dataset_type Name indicating which type of synthetic Visium data needs to be generated. \cr \cr
#' "real": for each region (defined in scRNAseq data - region_var - and thus 'real'), cell type frequencies are based on the number of cells in the corresponding region in the input scRNAseq data \cr \cr
#' "real_top1": for each region (defined in scRNAseq data and thus 'real'), we will only consider cell types that are most frequently present in that 1 region ("top 1"). Cell type frequencies are based on the number of cells in the corresponding region in the input scRNAseq data \cr \cr
#' "real_top1_uniform": for each region (defined in scRNAseq data - region_var - and thus 'real'), we will only consider cell types that are most frequently present in that 1 region ("top 1"). Cell type frequencies are equal for each cell type assigned to a region ("uniform") \cr \cr
#' "real_top2_overlap": for each region (defined in scRNAseq data - region_var - and thus 'real'), we will only consider cell types that if this particular region is among the 2 regions in which they are most abundant ("top 2"). As a consequence, there will be some "overlap" in cell type presence between two regions.  Cell type frequencies are based on the number of cells in the corresponding region in the input scRNAseq data \cr \cr
#' "real_top2_overlap_uniform": for each region (defined in scRNAseq data - region_var - and thus 'real'), we will only consider cell types that if this particular region is among the 2 regions in which they are most abundant ("top 2"). As a consequence, there will be some "overlap" in cell type presence between two regions. Cell type frequencies are equal for each cell type assigned to a region ("uniform") \cr \cr
#' "real_missing_celltypes_visium": the same as the 'real' dataset type, but now some random input scRNAseq cell types will be left out, and thus not used to generate synthetic data. \cr \cr
#' "artificial_uniform_distinct": celltype-prior region assignments are entirely artificial and thus not based on possibly present regional information in the input scRNAseq data. Each artifical prior region has a distinct set of cell types here, and cell type frequencies are the same for each cell type in each region. \cr \cr
#' "artificial_diverse_distinct": celltype-prior region assignments are entirely artificial and thus not based on possibly present regional information in the input scRNAseq data. Each artifical prior region has a distinct set of cell types here, and cell type frequencies are different for each cell type in each region.  \cr \cr
#' "artificial_uniform_overlap": celltype-prior region assignments are entirely artificial and thus not based on possibly present regional information in the input scRNAseq data. Each artifical prior region could have some overlap in cell type contribution with other regions, and cell type frequencies are the same for each cell type in each region.  \cr \cr
#' "artificial_diverse_overlap": celltype-prior region assignments are entirely artificial and thus not based on possibly present regional information in the input scRNAseq data. Each artifical prior region could have some overlap in cell type contribution with other regions, and cell type frequencies are different for each cell type in each region.  \cr \cr
#' "artificial_dominant_celltype_diverse": celltype-prior region assignments are entirely artificial and thus not based on possibly present regional information in the input scRNAseq data. Each artifical prior region could have some overlap in cell type contribution with other regions, and cell type frequencies are different for each cell type in each region. There is one dominant celltype that is present in each region, and is much more abundant than other cell types (5-15x more). \cr \cr
#' "artificial_partially_dominant_celltype_diverse: celltype-prior region assignments are entirely artificial and thus not based on possibly present regional information in the input scRNAseq data. Each artifical prior region could have some overlap in cell type contribution with other regions, and cell type frequencies are different for each cell type in each region. There is one dominant celltype that is present in all but one regions ("partially dominant"), and is much more abundant than other cell types, except for one region where it is equally abundant. \cr \cr
#' "artificial_missing_celltypes_visium": the same as the 'artificial_diverse_overlap' dataset type, but now some random input scRNAseq cell types will be left out, and thus not used to generate synthetic data. \cr \cr
#' "artificial_dominant_rare_celltype_diverse": the same as 'artificial_dominant_celltype_diverse', except that the dominant cell type that is present in all regions is much less abundant than other cell types. \cr \cr
#' "artificial_regional_rare_celltype_diverse": the same as 'artificial_dominant_rare_celltype_diverse', except that the rare cell type is now only present in one region instead of all regions.\cr \cr
#' "artificial_diverse_distinct_missing_celltype_sc": the same as 'artificial_diverse_distinct' except that one random cell type will be removed from the reference scRNAseq dataset at time of integration and evaluation. This to resemble a case where integration is done with an incomplete reference. \cr \cr
#' "artificial_diverse_overlap_missing_celltype_sc": the same as 'artificial_diverse_overlap' except that one random cell type will be removed from the reference scRNAseq dataset at time of integration and evaluation. This to resemble a case where integration is done with an incomplete reference. \cr \cr
#' @return A list with three sublists: \cr
#' region_assignments: gives the number of spots and cell type frequencies that will be used for each prior region \cr
#' dataset_properties: a tibble describing the dataset properties (based on the dataset_type argument) \cr
#' gold_standard_priorregion: a tibble indicating the presence and frequency of each scRNAseq-derived cell type in each prior region
#'
#' @import Seurat
#' @import dplyr
#' @import tibble
#' @import purrr
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' region_assignment = make_region_celltype_assignment(seurat_obj = seurat_obj, clust_var = "subclass", n_regions = 5, dataset_type = "artificial_missing_celltypes_visium")
#' }
#'
#' @export
#'
make_region_celltype_assignment = function(seurat_obj, clust_var, n_regions, dataset_type, region_var = NULL, dataset_id = "1", n_spots_min = 50, n_spots_max = 500){

  requireNamespace("dplyr")
  requireNamespace("Seurat")

  if(class(seurat_obj) != "Seurat"){
    stop("'seurat_obj' argument should be a Seurat Object")
  }
  if(!is.character(clust_var)){
    stop("clust_var argument should be a character vector with one element")
  }
  if(!clust_var %in% colnames(seurat_obj@meta.data)){
    stop("clust_var argument should be a valid column name of the seurat_obj meta.data")
  }
  if(!is.null(n_regions)){
    if(!is.numeric(n_regions)){
      stop("n_regions argument should be a numeric vector with one element or NULL if the region_var argument is defined")
    }
  }


  if(!dataset_type %in% c("real",
                            "real_top1",
                            "real_top1_uniform",
                            "real_top2_overlap",
                            "real_top2_overlap_uniform",
                            "real_missing_celltypes_visium",
                            "artificial_uniform_distinct",
                            "artificial_diverse_distinct",
                            "artificial_uniform_overlap",
                            "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse",
                            "artificial_partially_dominant_celltype_diverse",
                            "artificial_missing_celltypes_visium",
                          "artificial_dominant_rare_celltype_diverse",
                          "artificial_regional_rare_celltype_diverse",
                          "artificial_diverse_distinct_missing_celltype_sc",
                          "artificial_diverse_overlap_missing_celltype_sc")
  ){
    stop("the dataset_type argument is not valid. Check the documentation to see the different possibilities.")
  }
  if(!is.null(region_var)){
    if(!is.character(region_var)){
      stop("region_var argument should be NULL or a character vector with one element")
    }
    if(!region_var %in% colnames(seurat_obj@meta.data)){
      stop("region_var argument should be NULL or a valid column name of the seurat_obj meta.data")
    }
  } else {
    if(!dataset_type %in% c("artificial_uniform_distinct",
                            "artificial_diverse_distinct",
                            "artificial_uniform_overlap",
                            "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse",
                            "artificial_partially_dominant_celltype_diverse",
                            "artificial_missing_celltypes_visium",
                            "artificial_dominant_rare_celltype_diverse",
                            "artificial_regional_rare_celltype_diverse",
                            "artificial_diverse_distinct_missing_celltype_sc",
                            "artificial_diverse_overlap_missing_celltype_sc")
    ){
      stop("the dataset_type argument is not valid. If you dont provide a 'region_var', then you cannot use a dataset type that will use this region_var to define cell type frequencies. Check the documentation to see the different possibilities.")
    }
  }
  if(!is.character(dataset_id)){
    stop("dataset_id argument should be a character vector with one element")
  }
  if(!is.numeric(n_spots_min) | length(n_spots_min) != 1){
    stop("n_spots_min argument should be a numeric vector with one element")
  }
  if(!is.numeric(n_spots_max) | length(n_spots_max) != 1){
    stop("n_spots_max argument should be a numeric vector with one element")
  }

  seurat_obj$seurat_clusters_oi = droplevels(seurat_obj@meta.data[,clust_var] %>% as.factor())
  celltypes = unique(seurat_obj$seurat_clusters_oi)

  if(!is.null(region_var)){
    seurat_obj$prior_regions = droplevels(seurat_obj@meta.data[,region_var] %>% as.factor())
    regions = unique(seurat_obj$prior_regions)
    n_regions = length(regions)
  } else {
    regions = paste0("priorregion",seq(n_regions))
  }

  if(length(regions) >= length(celltypes)){
    stop("number of regions should be smaller than number of different cell types!")
  }

  # Generate the cell type frequencies per dataset type + return the properties of that dataset type
  if(dataset_type == "real"){ # for every region, we will dplyr::filter out the cells that belong to this region

    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      metadata_region_oi = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi) %>% dplyr::filter(prior_regions == region_oi)
      n_cells_region_oi = metadata_region_oi %>% dplyr::pull(cell_id) %>% unique() %>% length()

      present_celltypes_numbers = metadata_region_oi %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::count()
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dplyr::setdiff(celltypes, present_celltypes_numbers$seurat_clusters_oi %>% unique), n = 0)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/n_cells_region_oi)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)
    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "real", real_artificial = "real", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }

  if (dataset_type == "real_top1") { # put every cell type just with one region, then for every region, we will dplyr::filter out the filtered cells that belong to this region  - but keep the real proportions
    metadata = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi)

    filtered_region_celltype = metadata %>% dplyr::group_by(seurat_clusters_oi, prior_regions) %>% dplyr::count() %>%
      dplyr::ungroup()

    filtered_region_celltype = filtered_region_celltype %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::top_n(1, n) %>% dplyr::ungroup() %>% dplyr::select(-n)
    regions = filtered_region_celltype %>% dplyr::pull(prior_regions) %>% unique()
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% dplyr::filter(prior_regions == region_oi) %>% dplyr::pull(seurat_clusters_oi) %>% unique()
      metadata_region_oi = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi) %>% dplyr::filter(prior_regions == region_oi & seurat_clusters_oi %in% celltypes_oi)
      n_cells_region_oi = metadata_region_oi %>% dplyr::pull(cell_id) %>% unique() %>% length()

      present_celltypes_numbers = metadata_region_oi %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::count()
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dplyr::setdiff(celltypes, present_celltypes_numbers$seurat_clusters_oi %>% unique), n = 0)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/n_cells_region_oi)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)
    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "real_top1", real_artificial = "real", uniform_diverse = "diverse", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }


  if(dataset_type == "real_top1_uniform"){ #  # put every cell type just with one region, then for every region, we will dplyr::filter out the filtered cells that belong to this region  - make the proportions uniform
    metadata = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi)

    filtered_region_celltype = metadata %>% dplyr::group_by(seurat_clusters_oi, prior_regions) %>% dplyr::count() %>%
      dplyr::ungroup()

    filtered_region_celltype = filtered_region_celltype %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::top_n(1, n) %>% dplyr::ungroup() %>% dplyr::select(-n)
    regions = filtered_region_celltype %>% dplyr::pull(prior_regions) %>% unique()
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% dplyr::filter(prior_regions == region_oi) %>% dplyr::pull(seurat_clusters_oi) %>% unique()
      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = 1)
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/length(celltypes_oi))

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)
    })

    names(celltype_frequencies_list) = regions
    properties_df = tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "real_top1_uniform", real_artificial = "real", uniform_diverse = "uniform", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }
  if(dataset_type == "real_top2_overlap"){ # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types is according to real data
    metadata = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi)

    filtered_region_celltype = metadata %>% dplyr::group_by(seurat_clusters_oi, prior_regions) %>% dplyr::count() %>%
      dplyr::ungroup()

    filtered_region_celltype = filtered_region_celltype %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::top_n(2, n) %>% dplyr::ungroup() %>% dplyr::select(-n)
    regions = filtered_region_celltype %>% dplyr::pull(prior_regions) %>% unique()
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% dplyr::filter(prior_regions == region_oi) %>% dplyr::pull(seurat_clusters_oi) %>% unique()
      metadata_region_oi = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi) %>% dplyr::filter(prior_regions == region_oi & seurat_clusters_oi %in% celltypes_oi)
      n_cells_region_oi = metadata_region_oi %>% dplyr::pull(cell_id) %>% unique() %>% length()

      present_celltypes_numbers = metadata_region_oi %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::count()
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dplyr::setdiff(celltypes, present_celltypes_numbers$seurat_clusters_oi %>% unique), n = 0)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/n_cells_region_oi)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)
    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "real_top2_overlap", real_artificial = "real", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }

  if(dataset_type == "real_top2_overlap_uniform"){ #  # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types is uniform
    metadata = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi)

    filtered_region_celltype = metadata %>% dplyr::group_by(seurat_clusters_oi, prior_regions) %>% dplyr::count() %>%
      dplyr::ungroup()

    filtered_region_celltype = filtered_region_celltype %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::top_n(2, n) %>% dplyr::ungroup() %>% dplyr::select(-n)
    regions = filtered_region_celltype %>% dplyr::pull(prior_regions) %>% unique()
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% dplyr::filter(prior_regions == region_oi) %>% dplyr::pull(seurat_clusters_oi) %>% unique()
      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = 1)
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/length(celltypes_oi))

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)
    })

    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "real_top2_overlap_uniform", real_artificial = "real", uniform_diverse = "uniform", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }

  if (dataset_type == "real_missing_celltypes_visium"){ # take 4 missing cell types from scRNAseq data that are not used in Visium
    metadata = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi)

    filtered_region_celltype = metadata %>% dplyr::group_by(seurat_clusters_oi, prior_regions) %>% dplyr::count() %>%
      dplyr::ungroup()

    filtered_region_celltype = filtered_region_celltype %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::top_n(2, n) %>% dplyr::ungroup() %>% dplyr::select(-n)
    regions = filtered_region_celltype %>% dplyr::pull(prior_regions) %>% unique()

    n_leave_out_pool = min(length(celltypes) - 1,4)
    leave_out_celltypes = celltypes[sample(seq(length(celltypes)),n_leave_out_pool)]
    #print("cell types that will be left out: ")
    #print(leave_out_celltypes)
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% dplyr::filter(prior_regions == region_oi) %>% dplyr::pull(seurat_clusters_oi) %>% unique() %>%dplyr::setdiff(leave_out_celltypes)
      metadata_region_oi = seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, prior_regions, seurat_clusters_oi) %>% dplyr::filter(prior_regions == region_oi & seurat_clusters_oi %in% celltypes_oi)
      n_cells_region_oi = metadata_region_oi %>% dplyr::pull(cell_id) %>% unique() %>% length()

      present_celltypes_numbers = metadata_region_oi %>% dplyr::group_by(seurat_clusters_oi) %>% dplyr::count()
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dplyr::setdiff(celltypes, present_celltypes_numbers$seurat_clusters_oi %>% unique), n = 0)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/n_cells_region_oi)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)
      return(celltype_frequencies)
    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "real_missing_celltypes_visium", real_artificial = "real", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = TRUE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }

  if(dataset_type == "artificial_uniform_distinct"){ #   every region consists of different cell types, frequency uniform

    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)


    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]
      other_celltypes = setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = 1)
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/length(celltypes_oi))

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_uniform_distinct", real_artificial = "artificial", uniform_diverse = "uniform", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }
  if(dataset_type == "artificial_diverse_distinct"){ #  every region consists of different cell types, frequency random
    # celltype_region_splits =  split( celltypes, sample(n_regions, length(celltypes), repl = TRUE) )
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]

      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      sum_present = sum(present_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_diverse_distinct", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }
  if(dataset_type == "artificial_uniform_overlap"){ # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types will be uniform
    n_overlap_pool = min(length(celltypes),4)
    overlapping_celltypes = celltypes[sample(seq(length(celltypes)),n_overlap_pool)] %>% as.character() # and out of these 4, 1 will be chosen per region
    # celltype_region_splits =  split( celltypes, sample(n_regions, length(celltypes), repl = TRUE) )
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]
      celltypes_oi = union(celltypes_oi, sample(overlapping_celltypes,1))
      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = 1)
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/length(celltypes_oi))

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_uniform_overlap", real_artificial = "artificial", uniform_diverse = "uniform", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }
  if(dataset_type == "artificial_diverse_overlap"){ # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types will be random
    n_overlap_pool = min(length(celltypes),4)
    overlapping_celltypes = celltypes[sample(seq(length(celltypes)),n_overlap_pool)] %>% as.character() # and out of these 4, 1 will be chosen per region
    # celltype_region_splits =  split( celltypes, sample(n_regions, length(celltypes), repl = TRUE) )
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]
      celltypes_oi = union(celltypes_oi, sample(overlapping_celltypes,1))

      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      sum_present = sum(present_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_diverse_overlap", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }
  if (dataset_type == "artificial_dominant_celltype_diverse"){ # one cell type has significant higher frequency over all regions - difference in regions determined by other cell types. Some cell types in more than one region.

    dominant_celltype = sample(celltypes,1)
    #print("dominant cell type is")
    #print(dominant_celltype)

    # celltype_region_splits =  split( celltypes, sample(n_regions, length(celltypes), repl = TRUE) )
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]] %>%dplyr::setdiff(dominant_celltype)
      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi) %>%dplyr::setdiff(dominant_celltype)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)
      dominant_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dominant_celltype, n = sample(seq(75,100), 1))

      sum_present = sum(present_celltypes_numbers$n) %>% sum(dominant_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers) %>% dplyr::bind_rows(dominant_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)

    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_dominant_celltype_diverse", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = TRUE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }
  if (dataset_type == "artificial_partially_dominant_celltype_diverse"){ # one cell type has significant higher frequency over all regions - difference in regions determined by other cell types. Some cell types in more than one region.

    dominant_celltype = sample(celltypes,1)
    #print("dominant cell type is")
    #print(dominant_celltype)

    # celltype_region_splits =  split( celltypes, sample(n_regions, length(celltypes), repl = TRUE) )
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]] %>%dplyr::setdiff(dominant_celltype)
      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi) %>%dplyr::setdiff(dominant_celltype)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)
      if(i == length(regions)){
        dominant_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dominant_celltype, n = 0)
      } else if (i == length(regions) -1){
        dominant_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dominant_celltype, n = sample(seq(1,10), 1))
      } else {
        dominant_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dominant_celltype, n = sample(seq(75,100), 1))
      }

      sum_present = sum(present_celltypes_numbers$n) %>% sum(dominant_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers) %>% dplyr::bind_rows(dominant_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)

    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_partially_dominant_celltype_diverse", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = TRUE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)

  }


  if (dataset_type == "artificial_missing_celltypes_visium"){ # take 4 missing cell types from scRNAseq data that are not used in Visium
    n_overlap_pool = min(length(celltypes),4)
    overlapping_celltypes = celltypes[sample(seq(length(celltypes)),n_overlap_pool)] %>% as.character() # and out of these 4, 1 will be chosen per region

    n_leave_out_pool = min(length(celltypes) - 1,4)
    leave_out_celltypes = celltypes[sample(seq(length(celltypes)),n_leave_out_pool)] %>% as.character()
    #print("cell types that will be left out: ")
    #print(leave_out_celltypes)

    overlapping_celltypes = overlapping_celltypes %>% dplyr::setdiff(leave_out_celltypes)

    # celltype_region_splits =  split( celltypes, sample(n_regions, length(celltypes), repl = TRUE) )
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]
      celltypes_oi = union(celltypes_oi, sample(overlapping_celltypes,1))
      celltypes_oi = dplyr::setdiff(celltypes_oi, leave_out_celltypes)
      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      sum_present = sum(present_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_missing_celltypes_visium", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = TRUE, rare_celltype = FALSE, missing_celltype_sc = NA)
  }

  if (dataset_type == "artificial_diverse_distinct_missing_celltype_sc"){ #  every region consists of different cell types, frequency random

    missing_celltype_sc = sample(celltypes,1)

    # celltype_region_splits =  split( celltypes, sample(n_regions, length(celltypes), repl = TRUE) )
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]

      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      sum_present = sum(present_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_diverse_distinct_missing_celltype_sc", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = missing_celltype_sc)

  }
  if (dataset_type == "artificial_diverse_overlap_missing_celltype_sc"){ # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types will be random

    missing_celltype_sc = sample(celltypes,1)

    n_overlap_pool = min(length(celltypes),4)
    overlapping_celltypes = celltypes[sample(seq(length(celltypes)),n_overlap_pool)] %>% as.character() # and out of these 4, 1 will be chosen per region
    # celltype_region_splits =  split( celltypes, sample(n_regions, length(celltypes), repl = TRUE) )
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(celltypes, celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]
      celltypes_oi = union(celltypes_oi, sample(overlapping_celltypes,1))

      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)

      sum_present = sum(present_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)

    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_diverse_overlap_missing_celltype_sc", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = missing_celltype_sc)

  }
  if (dataset_type == "artificial_dominant_rare_celltype_diverse"){ # one cell type that is present over all regions, but in low frequency - difference in regions determined by other cell types. Some cell types in more than one region.

    dominant_celltype = sample(celltypes,1)
    #print("dominant cell type is")
    #print(dominant_celltype)

    # check that the rare cell type is never alone in one region
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) -1 - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(setdiff(celltypes,dominant_celltype), celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]] %>%dplyr::setdiff(dominant_celltype)
      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi) %>%dplyr::setdiff(dominant_celltype)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(seq(75,125),length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)
      dominant_celltypes_numbers = tibble::tibble(seurat_clusters_oi = dominant_celltype, n = 10)

      sum_present = sum(present_celltypes_numbers$n) %>% sum(dominant_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers) %>% dplyr::bind_rows(dominant_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)

    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_dominant_rare_celltype_diverse", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = TRUE, missing_celltype = FALSE, rare_celltype = TRUE, missing_celltype_sc = NA)

  }
  if (dataset_type == "artificial_regional_rare_celltype_diverse"){ # one cell type that is present over all regions, but in low frequency - difference in regions determined by other cell types. Some cell types in more than one region.

    rare_celltype = sample(celltypes,1)
    region_of_rare = sample(regions,1)
    #print("rare cell type is")
    #print(rare_celltype)

    # check that the rare cell type is never alone in one region
    celltype_region_samplings = c(sample(seq(n_regions)), sample(seq(n_regions), size = length(celltypes) -1 - length(seq(n_regions)), replace = TRUE))
    celltype_region_splits =  split(setdiff(celltypes,rare_celltype), celltype_region_samplings)

    celltype_frequencies_list = seq(n_regions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]] %>%dplyr::setdiff(rare_celltype)
      other_celltypes = dplyr::setdiff(celltypes, celltypes_oi) %>%dplyr::setdiff(rare_celltype)

      present_celltypes_numbers = tibble::tibble(seurat_clusters_oi = celltypes_oi, n = sample(seq(75,125),length(celltypes_oi)))
      absent_celltypes_numbers = tibble::tibble(seurat_clusters_oi = other_celltypes, n = 0)
      if(region_oi == region_of_rare){
        rare_celltypes_numbers = tibble::tibble(seurat_clusters_oi = rare_celltype, n = 10)
      } else {
        rare_celltypes_numbers = tibble::tibble(seurat_clusters_oi = rare_celltype, n = 0)
      }

      sum_present = sum(present_celltypes_numbers$n) %>% sum(rare_celltypes_numbers$n)

      celltype_numbers = dplyr::bind_rows(present_celltypes_numbers, absent_celltypes_numbers) %>% dplyr::bind_rows(rare_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% dplyr::mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% dplyr::pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% dplyr::pull(seurat_clusters_oi)

      return(celltype_frequencies)

    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble::tibble(dataset_id = paste0(dataset_type, dataset_id),dataset_type = "artificial_regional_rare_celltype_diverse", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = TRUE, missing_celltype_sc = NA)

  }

  # for in next function, we will need to have the variable in which region is stored if generating 'real synthetic data'
  properties_df = properties_df %>% dplyr::inner_join(tibble::tibble(real_artificial = c("real","artificial"), real_region_var = c(region_var, NA)), by = "real_artificial")
  # print(properties_df)
  # Define the number of spots that needs to be generated per region (random number between n_spots_min and n_spots_max)
  n_spots_list = regions %>% lapply(function(region_oi){
    return(sample(seq(n_spots_min,n_spots_max),1))
  })
  names(n_spots_list) = regions

  celltype_frequencies_list = celltype_frequencies_list[names(n_spots_list)]
  region_assignments = purrr::map2(n_spots_list, celltype_frequencies_list, function(x,y){return(list(n_spots = x, celltype_freq = y))})
  gs = purrr::map2(names(celltype_frequencies_list),celltype_frequencies_list, function(region, celltype_freq){
    tibble::tibble(prior_region = region, celltype = names(celltype_freq), freq = celltype_freq, present = celltype_freq > 0)
  }) %>% dplyr::bind_rows()

  return(list(region_assignments = region_assignments, dataset_properties = properties_df, gold_standard_priorregion = gs))

}
#' @title Generate synthetic Visium/ST data from input scRNAseq data based on predefined cell type frequencies and number of spots.
#'
#' @description \code{region_assignment_to_syn_data} Generate synthetic Visium/ST data based on predefined cell type frequencies and number of spots.
#' To create synthetic Visium spots, we will generate an artificial mixture of single cells (based on the assumption that a spot is a mixture of n cells (n between 2 and 10))
#' Cells from the input scRNAseq data will be sampled based on the provided cell type frequencies per region (such that each region has a different cellular composition), and their counts will be summed.
#' After summation of the counts, the counts will be downsampled to get a total number of counts per spot that is typical for a Visium dataset. This target number of counts will be determined by sampling from a Normal distribution with mean and SD provided by the user.
#' @usage
#' region_assignment_to_syn_data(region_assignment_list, seurat_obj, clust_var, visium_mean = 20000, visium_sd = 7000, n_cells_min = 2, n_cells_max = 10, add_mock_region = FALSE)
#'
#' @inheritParams make_region_celltype_assignment
#' @param region_assignment_list Output of the function `make_region_celltype_assignment`
#' @param visium_mean The mean of the normal distribution that will be used to pick the target number of counts for downsample.
#' @param visium_sd The standard deviation of the normal distribution that will be used to pick the target number of counts for downsample.
#' @param n_cells_min Integer indicating the minimum number of cells sampled per spot. Default: 2
#' @param n_cells_max  Integer indicating the maximum number of cells sampled per spot. Default: 10
#' @param add_mock_region If you want to make a "mock region" in addition to the other regions, set this parameter to TRUE. Mock region is recommended when the synthetic data will be used to evaluate spatial/region annotation of cells, not when evaluating deconvolution tools. Default: FALSE. The mock region is generated similar as the real regions, but the input cell type frequencies are the same for each cell type, and after generation of the counts, the gene names are shuffled such that cellular identities in this mockregion are lost.
#'
#' @return list with five sublists: \cr
#' counts: count matrix of the synthetic visium data \cr
#' spot_composition: the cell type composition of each spot \cr
#' relative_spot_composition: the relative cell type composition of each spot \cr
#' gold_standard_priorregion: indication of which cell types belong to which prior region \cr
#' dataset_properties: indication of the properties of the dataset (see `make_region_celltype_assignment`)
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' region_assignment = make_region_celltype_assignment(seurat_obj = seurat_obj, clust_var = "subclass", n_regions = 5, dataset_type = "artificial_missing_celltypes_visium")
#' synthetic_visium_data = region_assignment_to_syn_data(region_assignment_list = region_assignment, seurat_obj = seurat_obj, clust_var = "subclass")
#' }
#'
#' @export
#'
region_assignment_to_syn_data = function(region_assignment_list, seurat_obj, clust_var, visium_mean = 20000, visium_sd = 7000,
                                         n_cells_min = 2, n_cells_max = 10, add_mock_region = FALSE){

  requireNamespace("dplyr")
  requireNamespace("Seurat")

  if(class(seurat_obj) != "Seurat"){
    stop("'seurat_obj' argument should be a Seurat Object")
  }
  if(!is.list(region_assignment_list)){
    stop("region_assignment_list should be a list")
  }
  if(! "region_assignments" %in% names(region_assignment_list)){
    stop("region_assignment_list should be a list with 'region_assignments' as sublist")
  }
  if(! "dataset_properties" %in% names(region_assignment_list)){
    stop("region_assignment_list should be a list with 'dataset_properties' as sublist")
  }
  if(! "gold_standard_priorregion" %in% names(region_assignment_list)){
    stop("region_assignment_list should be a list with 'gold_standard_priorregion' as sublist")
  }
  if(!is.character(clust_var)){
    stop("clust_var argument should be a character vector with one element")
  }
  if(!clust_var %in% colnames(seurat_obj@meta.data)){
    stop("clust_var argument should be a valid column name of the seurat_obj meta.data")
  }
  if(!is.numeric(visium_mean) | length(visium_mean) != 1){
    stop("visium_mean argument should be a numeric vector with one element")
  }
  if(!is.numeric(visium_sd) | length(visium_sd) != 1){
    stop("visium_sd argument should be a numeric vector with one element")
  }

  region_assignments = region_assignment_list[["region_assignments"]]
  dataset_properties = region_assignment_list[["dataset_properties"]]
  real_dataset = dataset_properties$real_artificial == "real"
  gold_standard_priorregion = region_assignment_list[["gold_standard_priorregion"]]

  seurat_obj$seurat_clusters_oi = droplevels(seurat_obj@meta.data[,clust_var] %>% as.factor())
  if(real_dataset == TRUE){
    region_var = dataset_properties$real_region_var
    seurat_obj$prior_regions = droplevels(seurat_obj@meta.data[,region_var] %>% as.factor())
  }

  regions = names(region_assignments)
  synthetic_regions = regions %>% lapply(generate_spots, region_assignments, seurat_obj, visium_mean, visium_sd, 
                                         n_cells_min, n_cells_max, real_dataset)

  if(add_mock_region == TRUE){   # possible to add the mock region here in the synthetic_regions thing
    mock_regions = "mockregion"
    n_spots_mock = region_assignments %>% purrr::map("n_spots") %>% unlist() %>% min()
    celltypes = region_assignments %>% purrr::map("celltype_freq") %>% .[[1]] %>% names()
    celltype_freq_mock = rep(1/length(celltypes), times = length(celltypes))
    names(celltype_freq_mock) = celltypes
    region_assignments_mock_region = list(n_spots = n_spots_mock,
                                          celltype_freq = celltype_freq_mock)
    region_assignments_mock = list(region_assignments_mock_region)
    names(region_assignments_mock) = mock_regions

    synthetic_regions_mock = mock_regions %>% lapply(generate_spots, region_assignments_mock, seurat_obj, visium_mean, visium_sd,
                                                     n_cells_min, n_cells_max, real_dataset = FALSE)
    # synthetic_regions_mock[[1]]$spot_composition %>% select(-name, -region) %>% apply(2,sum)

    # now shuffle the gene names at random
    mock_genes = rownames(synthetic_regions_mock[[1]]$counts)
    shuffled_mock_genes = sample(mock_genes, size = length(mock_genes), replace = F)
    rownames(synthetic_regions_mock[[1]]$counts) = shuffled_mock_genes
    # but keep the original order
    synthetic_regions_mock[[1]]$counts = synthetic_regions_mock[[1]]$counts %>% .[mock_genes,]

    mock_spot_composition = synthetic_regions_mock[[1]]$spot_composition
    mock_spot_composition[mock_spot_composition > 0] = 0
    mock_spot_composition$name = synthetic_regions_mock[[1]]$spot_composition$name
    mock_spot_composition$region =  synthetic_regions_mock[[1]]$spot_composition$region
    synthetic_regions_mock[[1]]$spot_composition = mock_spot_composition

    synthetic_regions = c(synthetic_regions, synthetic_regions_mock)

    gold_standard_priorregion_mock = tibble(prior_region = mock_regions, celltype = celltypes, freq = 0, present = FALSE)
    gold_standard_priorregion = gold_standard_priorregion %>% bind_rows(gold_standard_priorregion_mock)
  }

  syn_count_data = synthetic_regions %>% purrr::map("counts") %>% do.call("cbind",.)
  syn_spot_composition = synthetic_regions %>% purrr::map("spot_composition") %>% dplyr::bind_rows()
  syn_spot_composition_relative = syn_spot_composition %>% .[,-which(colnames(.) %in% c("name","region"))] %>% apply(1, function(x){x/sum(x)}) %>% t() %>% cbind(.,syn_spot_composition[,c("name","region")])
  return(list(counts = syn_count_data, spot_composition = syn_spot_composition, relative_spot_composition = syn_spot_composition_relative, gold_standard_priorregion = gold_standard_priorregion, dataset_properties = dataset_properties))

}
#' @title Generate the synthetic data of a spot from a prior region of interest.
#'
#' @description \code{generate_spots} Generate the synthetic data of a spot from a prior region of interest.
#' To create a synthetic Visium spot, we will generate an artificial mixture of single cells (based on the assumption that a spot is a mixture of n cells)
#' Cells from the input scRNAseq data will be sampled based on the provided cell type frequencies of the region of interest, and their counts will be summed.
#' After summation of the counts, the counts will be downsampled to get a total number of counts per spot that is typical for a Visium dataset. This target number of counts will be determined by sampling from a Normal distribution with mean and SD provided by the user.
#' @usage
#' generate_spots(region_oi, region_assignments, seurat_obj, visium_mean = 20000, visium_sd = 7000, n_cells_min = 2, n_cells_max = 10, real_dataset = FALSE)
#'
#' @inheritParams region_assignment_to_syn_data
#' @inheritParams make_region_celltype_assignment
#' @param region_oi Name of the prior region of interest from which the spot data should be generated.
#' @param region_assignments List indicating the cell type frequencies for each prior region. This could be the `region_assignments` sublist of the output of `make_region_celltype_assignment`
#' @param real_dataset Logical vector indicating whether the synthetic dataset is based on sampling from 'real' regions. See the documentation of the `make_region_celltype_assignment` function for more information about the distinction between 'real' and 'artificial' synthetic datasets. Default: FALSE.
#'
#' @return A list with two sublists: the counts and the cell type composition of the spot.
#'
#' @import Seurat
#' @import dplyr
#' @import dtplyr
#' @import tidyr
#' @importFrom Matrix Matrix
#' @importFrom  DropletUtils downsampleMatrix
#' @importFrom stats rnorm
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(Seurat)
#' region_assignment_list = make_region_celltype_assignment(seurat_obj = seurat_obj, clust_var = "subclass", n_regions = 5, dataset_type = "artificial_missing_celltypes_visium")
#' region_assignments = region_assignment_list[["region_assignments"]]
#' region_oi = names(region_assignments) %>% head(1)
#' spot = generate_spots(region_oi = region_oi, region_assignments = region_assignments, seurat_obj = seurat_obj)
#' }
#'
#' @export
#'
generate_spots = function(region_oi, region_assignments, seurat_obj, visium_mean = 20000, visium_sd = 7000,
                          n_cells_min = 2, n_cells_max = 10, real_dataset = FALSE){

  requireNamespace("dplyr")
  requireNamespace("dtplyr") # use dplyr commands with data.table speed
  requireNamespace("Seurat")

  if(class(seurat_obj) != "Seurat"){
    stop("'seurat_obj' argument should be a Seurat Object")
  }
  if(!is.list(region_assignments)){
    stop("region_assignments should be a list")
  }
  if(!is.numeric(visium_mean) | length(visium_mean) != 1){
    stop("visium_mean argument should be a numeric vector with one element")
  }
  if(!is.numeric(visium_sd) | length(visium_sd) != 1){
    stop("visium_sd argument should be a numeric vector with one element")
  }
  if(!is.logical(real_dataset)){
    stop("real_dataset should be TRUE or FALSE")
  }
  if(!is.character(region_oi)){
    stop("region_oi should be a character vector")
  } else {
    if(! region_oi %in% names(region_assignments)){
      stop("region_oi should be a name of a prior region as defined in 'region_assignments")
    }
  }

  count_matrix = as.matrix(seurat_obj[["RNA"]]@counts)
  metadata = seurat_obj@meta.data
  metadata = metadata %>% tibble::rownames_to_column("cell_id") %>% tibble::as_tibble()

  info_oi = region_assignments[[region_oi]]
  n = info_oi$n_spots
  celltype_freq = info_oi$celltype_freq

  # #print(region_oi)
  # #print(celltype_freq)
  ds_spots = lapply(1:n, function(i){

    # select between n_cells_min and n_cells_max cells (default: 2, 10) randomly from the count matrix
    # add here something that uses the cell type frequencies!
    n_cells = sample(x = n_cells_min:n_cells_max, size = 1)

    celltypes_pool = sample(x = names(celltype_freq), size =  n_cells, prob = celltype_freq, replace = T)
    celltypes_pool_table = celltypes_pool %>% table()

    cell_pool = names(celltypes_pool_table) %>% lapply(function(celltype_oi, celltypes_pool_table, metadata, real_dataset){
      n_cells_oi = celltypes_pool_table[celltype_oi]
      # #print(n_cells_oi)
      if(real_dataset == FALSE){
        source_cells = metadata %>% dplyr::filter(seurat_clusters_oi == celltype_oi) %>% dplyr::pull(cell_id) #
      } else {
        source_cells = metadata %>% dplyr::filter(seurat_clusters_oi == celltype_oi & prior_regions == region_oi) %>% dplyr::pull(cell_id) # if dataset is real: pick only cells from possible region! (but then it should be present right...)
      }
      # #print(length(source_cells))
      picked_cells = sample(source_cells, size = n_cells_oi, replace = T)
      return(picked_cells)
    }, celltypes_pool_table, metadata, real_dataset) %>% unlist() %>% unique()

    # Create a name for the spot with the necessary info to deconvolute it
    pos = which(colnames(count_matrix) %in% cell_pool)
    tmp_ds = seurat_obj@meta.data[pos,] %>% dplyr::mutate(weight = 1)
    name_simp = paste(region_oi,'_spot_',i,sep='')

    # table with information of cell types per spots
    spot_ds = suppressMessages(tmp_ds %>%
                                 dplyr::select(seurat_clusters_oi, weight) %>%
                                 dplyr::mutate(seurat_clusters_oi = as.character(seurat_clusters_oi)) %>%
                                 # dplyr::mutate(seurat_clusters_oi = paste('clust_',seurat_clusters_oi,sep = '')) %>%
                                 dplyr::group_by(seurat_clusters_oi) %>%
                                 dplyr::summarise(sum_weights = sum(weight)) %>% #
                                 tidyr::pivot_wider(names_from = seurat_clusters_oi, values_from = sum_weights) %>%
                                 dplyr::mutate(name = name_simp, region = region_oi))

    # Generate synthetic spot

    ## Here we add up the counts of each cell
    syn_spot = rowSums(as.matrix(count_matrix[,cell_pool])); sum(syn_spot)
    names_genes = names(syn_spot)

    ## Downsample
    ### 25k is a bit above average 20k UMIs observed in spatial transcriptomics data then downsample to 20k: spotlight
    ### we do it as such: if counts higher than mean + sd: sample again: we

    if(sum(syn_spot) > (visium_mean + visium_sd)){
      prop_oi = 0
      while(prop_oi <= 0){
        prop_oi = floor(rnorm(1,visium_mean, visium_sd))/sum(syn_spot)
      }
      if(prop_oi > 1){
        prop_oi = 1
      }

      syn_spot_sparse = DropletUtils::downsampleMatrix(Matrix::Matrix(syn_spot,sparse = T), prop = prop_oi)
    } else {
      syn_spot_sparse = Matrix::Matrix(syn_spot,sparse = T)
    }

    rownames(syn_spot_sparse) = names_genes
    colnames(syn_spot_sparse) = name_simp

    # Return the transpose so that spots are on the rows and genes are on the columns
    return(list(syn_spot_sparse,spot_ds))
  })
  ds_syn_spots = map(ds_spots, 1) %>%
    base::Reduce(function(m1,m2) cbind(unlist(m1),unlist(m2)), .)

  # Generate dataframe of spot characteristic
  ds_spots_metadata = map(ds_spots, 2) %>%
    dplyr::bind_rows() %>%
    dtplyr::lazy_dt(.) %>%  # Convert to lazy data so it can be avaluated as a data.table
    # dplyr::mutate_all(~replace_na(.,0)) %>%
    data.frame()

  ds_spots_metadata[is.na(ds_spots_metadata)] = 0
  # change column order so that its progressive
  lev_mod = gsub("[\\+|\\ ]", ".", levels(seurat_obj$seurat_clusters_oi))
  all_cn = c(lev_mod,'name','region')
  all_cn = make.names(all_cn) # needed for when not all cell type names are suitable column names according to R
  if( sum(all_cn %in% colnames(ds_spots_metadata)) == (nlevels(seurat_obj$seurat_clusters_oi)+2) ){
    ds_spots_metadata = ds_spots_metadata[,all_cn]
  } else {
    missing_cols = all_cn[ which( ! all_cn %in% colnames(ds_spots_metadata) ) ]
    ds_spots_metadata[missing_cols] = 0
    ds_spots_metadata = ds_spots_metadata[,all_cn]
  }
  return(list(counts = ds_syn_spots, spot_composition = ds_spots_metadata))

}
#' @title Generate synthetic Visium/ST data from input scRNAseq data.
#'
#' @description \code{generate_synthetic_visium} This function generates synthetic Visium/ST data from input scRNAseq data. \cr
#' To create synthetic Visium spots, we will generate an artificial mixture of single cells (based on the assumption that a spot is a mixture of n cells \cr
#' Cells from the input scRNAseq data will be sampled based on the provided cell type frequencies per region (such that each region has a different cellular composition), and their counts will be summed. \cr
#' After summation of the counts, the counts will be downsampled to get a total number of counts per spot that is typical for a Visium dataset. This target number of counts will be determined by sampling from a Normal distribution with mean and SD provided by the user. \cr
#' You can generate different 'tissue types' by changing the way the cell type frequencies are defined. Read the argument description of `dataset_type` for seeing the different tissue type options. \cr
#' If your scRNAseq data object already provides information about the region/location a cell was sampled from, you could use this information as well to generate synthetic Visium data that will only pool cells together if they originate from the same region. \cr
#'
#' @usage
#' generate_synthetic_visium(seurat_obj, dataset_type, clust_var, n_regions, region_var = NULL, dataset_id = "1", n_spots_min = 50, n_spots_max = 500, visium_mean = 20000, visium_sd = 7000, n_cells_min = 2, n_cells_max = 10, add_mock_region = FALSE, sc_rnaseq_path = NA)
#'
#' @inheritParams make_region_celltype_assignment
#' @inheritParams region_assignment_to_syn_data
#' @param sc_rnaseq_path NA or indication of the path to read in the scRNAseq that should be used to integratie with this synthetic visium data. Default: NA
#'
#' @return list with six sublists: \cr
#' counts: count matrix of the synthetic visium data \cr
#' spot_composition: the cell type composition of each spot \cr
#' relative_spot_composition: the relative cell type composition of each spot \cr
#' gold_standard_priorregion: indication of which cell types belong to which prior region \cr
#' dataset_properties: indication of the properties of the dataset
#' sc_rnaseq_path: NA or indication of the path to read in the scRNAseq that should be used to integratie with this synthetic visium data
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_missing_celltypes_visium", clust_var = "subclass", n_regions = 5)
#' }
#'
#' @export
#'
generate_synthetic_visium = function(seurat_obj, dataset_type, clust_var, n_regions, region_var = NULL, dataset_id = "1", n_spots_min = 50, n_spots_max = 500, visium_mean = 20000, visium_sd = 7000, n_cells_min = 2, n_cells_max = 10, add_mock_region = FALSE, sc_rnaseq_path = NA){

  requireNamespace("dplyr")
  requireNamespace("Seurat")

  # input checks are implemented in the functions that are used here under the hood.

  region_assignment_list = make_region_celltype_assignment(seurat_obj = seurat_obj, clust_var = clust_var, n_regions = n_regions, dataset_type = dataset_type, region_var = region_var, dataset_id = dataset_id, n_spots_min = n_spots_min, n_spots_max = n_spots_max)
  synthetic_visium_data = region_assignment_to_syn_data(region_assignment_list = region_assignment_list, seurat_obj = seurat_obj, clust_var = clust_var, visium_mean = visium_mean, visium_sd = visium_sd,
                                                        n_cells_min = n_cells_min, n_cells_max = n_cells_max, add_mock_region = add_mock_region)
  synthetic_visium_data$sc_rnaseq_path = sc_rnaseq_path
  return(synthetic_visium_data)

}

#' @title Generate synthetic Visium/ST data from input scRNAseq data, using the cell type composition from the input as priors.
#'
#' @description \code{generate_synthetic_visium_lite} This function generates synthetic Visium/ST data from input scRNAseq data. \cr
#' In contrast to the original function, the cell type frequencies are based on the cell type composition of the input scRNAseq data. \cr
#' Because of this, there can only be 1 region. It is recommended to use single-nucleus RNA-seq data, as it has been shown to capture tissue compositions more accurately than scRNA-seq. \cr
#'
#' @usage
#' generate_synthetic_visium_lite(seurat_obj, clust_var, dataset_id = "1", n_spots = 500, visium_mean = 20000, visium_sd = 7000, n_cells_min = 2, n_cells_max = 10, sc_rnaseq_path = NA)
#'
#' @param seurat_obj The input scRNAseq data stored as Seurat object.
#' @param clust_var Name of the meta data column in which the cell type label of interest is defined.
#' @param dataset_id Character vector to define how you want to identify the dataset. The final identifier/name of the dataset will be the description of the dataset type followed by this character vector.
#' @param n_spots The number of spots to be generated. Default 500.
#' @param sc_rnaseq_path NA or indication of the path to read in the scRNAseq that should be used to integratie with this synthetic visium data. Default: NA
#' @param visium_mean The mean of the normal distribution that will be used to pick the target number of counts for downsample.
#' @param visium_sd The standard deviation of the normal distribution that will be used to pick the target number of counts for downsample.
#' @param n_cells_min Integer indicating the minimum number of cells sampled per spot. Default: 2
#' @param n_cells_max  Integer indicating the maximum number of cells sampled per spot. Default: 10
#' 
#' @return list with six sublists: \cr
#' counts: count matrix of the synthetic visium data \cr
#' spot_composition: the cell type composition of each spot \cr
#' relative_spot_composition: the relative cell type composition of each spot \cr
#' gold_standard_priorregion: indication of which cell types belong to which prior region \cr
#' dataset_properties: indication of the properties of the dataset
#' sc_rnaseq_path: NA or indication of the path to read in the scRNAseq that should be used to integratie with this synthetic visium data
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' synthetic_visium_data = generate_synthetic_visium_lite(seurat_obj = seurat_obj, clust_var = "subclass")
#' }
#'
#' @export
#'
generate_synthetic_visium_lite = function(seurat_obj, clust_var, dataset_id = "1", n_spots = 500, visium_mean = 20000, visium_sd = 7000,
                                          n_cells_min = 2, n_cells_max = 10, sc_rnaseq_path = NA){
  
  requireNamespace("dplyr")
  requireNamespace("Seurat")
  
  if(class(seurat_obj) != "Seurat"){
    stop("'seurat_obj' argument should be a Seurat Object")
  }
  if(!is.character(clust_var)){
    stop("clust_var argument should be a character vector with one element")
  }
  if(!clust_var %in% colnames(seurat_obj@meta.data)){
    stop("clust_var argument should be a valid column name of the seurat_obj meta.data")
  }
  if(!is.character(dataset_id)){
    stop("dataset_id argument should be a character vector with one element")
  }
  if(!is.numeric(n_spots) | length(n_spots) != 1){
    stop("n_spots argument should be a numeric vector with one element")
  }
  
  # Get cell type frequency based on cluster metadata
  seurat_obj$seurat_clusters_oi = droplevels(seurat_obj@meta.data[,clust_var] %>% as.factor())
  celltype_frequencies_list = list(priorregion1 = table(seurat_obj$seurat_clusters_oi)/ncol(seurat_obj))
  
  # Dataset property
  properties_df =  tibble::tibble(dataset_id = paste0("prior_from_data", dataset_id),dataset_type = "prior_from_data", real_artificial = "NA", uniform_diverse = "NA", distinct_overlap = "NA", dominant_celltype = FALSE, missing_celltype = FALSE, rare_celltype = FALSE, missing_celltype_sc = NA)
  
  # There is only 1 region, but create a list to be consistent with previous functions
  n_spots_list = list(n_spots) %>% setNames("priorregion1")
  region_assignments = purrr::map2(n_spots_list, celltype_frequencies_list, function(x,y){return(list(n_spots = x, celltype_freq = y))})
  gs = purrr::map2(names(celltype_frequencies_list),celltype_frequencies_list, function(region, celltype_freq){
    tibble::tibble(prior_region = region, celltype = names(celltype_freq), freq = celltype_freq, present = celltype_freq > 0)
  }) %>% dplyr::bind_rows()
  
  # Same workflow as original generate_synthetic_visium function
  region_assignment_list = list(region_assignments = region_assignments, dataset_properties = properties_df, gold_standard_priorregion = gs)
  synthetic_visium_data = region_assignment_to_syn_data(region_assignment_list = region_assignment_list, seurat_obj = seurat_obj, clust_var = clust_var, 
                                                        visium_mean = visium_mean, visium_sd = visium_sd, n_cells_min = n_cells_min, n_cells_max = n_cells_max)
  synthetic_visium_data$sc_rnaseq_path = sc_rnaseq_path
  
  return(synthetic_visium_data)
  
}
