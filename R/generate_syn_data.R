#' @title Get a Seurat object only containing cells of a region (or regions) of interest.
#'
#' @description \code{get_seuratObj_region_seurat} Get a Seurat object only containing cells of a region (or regions) of interest. This filtering is based on spatial annotation of cells as done by `predict_regions_of_cells`, and via the Seurat Label Transfer method.
#' @usage
#' get_seuratObj_region_seurat(region_oi, seurat_obj_scrnaseq, spatial_annotation_assay = "regions", label_transfer_cutoff = 0.25)
#'
#' @param seurat_obj Name of the assay of the Seurat object where the spatial region annotation predictions of cells were stored. Default: 'regions', as resulting from running `predict_regions_of_cells`
#' @param clust_var Cells with a spatial annotation score larger or equal than this cutoff will be kept, and thus considered as belonging to the region of interest. Default = 0.25.
#' @param nregions Cells with a spatial annotation score larger or equal than this cutoff will be kept, and thus considered as belonging to the region of interest. Default = 0.25.#' @param clust_var Cells with a spatial annotation score larger or equal than this cutoff will be kept, and thus considered as belonging to the region of interest. Default = 0.25.
#' @param dataset_type Cells with a spatial annotation score larger or equal than this cutoff will be kept, and thus considered as belonging to the region of interest. Default = 0.25.
#' @param region_var Cells with a spatial annotation score larger or equal than this cutoff will be kept, and thus considered as belonging to the region of interest. Default = 0.25.
#' @param dataset_name Cells with a spatial annotation score larger or equal than this cutoff will be kept, and thus considered as belonging to the region of interest. Default = 0.25.
#'
#' @return The Seurat object with the scRNAseq data (seurat_obj_scrnaseq) of cells belonging to the region(s) of interest.
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(Seurat)
#' seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
#' seurat_obj_scrnaseq_central = get_seuratObj_region_seurat(region_oi = "central", seurat_obj_scrnaseq, spatial_annotation_assay = "regions", label_transfer_cutoff = 0.25)
#' }
#'
#' @export
#'
make_region_celltype_assignment = function(seurat_obj, clust_var, nregions, dataset_type, region_var = NULL, dataset_name = "1"){

  requireNamespace("dplyr")
  requireNamespace("Seurat")

  if(class(seurat_obj_scrnaseq) != "Seurat")
    stop("'seurat_obj_scrnaseq' should be a Seurat Object")

  # dataset_types = c("real_sampling","real_sampling_top1","real_sampling_top1_uniform","real_sampling_top2_overlap","real_sampling_top2_overlap_uniform", "real_missing_celltypes_visium",
  #                   "artificial_sampling_uniform_distinct", "artificial_sampling_diverse_distinct", "artificial_sampling_uniform_overlap", "artificial_sampling_diverse_overlap", "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse")

  seurat_obj$seurat_clusters = droplevels(seurat_obj@meta.data[,clust_var] %>% as.factor())
  celltypes = unique(seurat_obj$seurat_clusters)

  if(!is.null(region_var)){
    seurat_obj$prior_regions = droplevels(seurat_obj@meta.data[,region_var] %>% as.factor())
    regions = unique(seurat_obj$prior_regions)
    nregions = length(regions)
  } else {
    regions = paste0("priorregion",seq(nregions))
  }

  if(length(regions) >= length(celltypes)){
    stop("number of regions should be smaller than number of different cell types!")
  }

  # Generate the cell type frequencies per dataset type + return the properties of that datast type
  if(dataset_type == "real_sampling"){ # for every region, we will filter out the cells that belong to this region

    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      metadata_region_oi = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters) %>% filter(prior_regions == region_oi)
      n_cells_region_oi = metadata_region_oi %>% pull(cell_id) %>% unique() %>% length()

      present_celltypes_numbers = metadata_region_oi %>% group_by(seurat_clusters) %>% dplyr::count()
      absent_celltypes_numbers = tibble(seurat_clusters = setdiff(celltypes, present_celltypes_numbers$seurat_clusters %>% unique), n = 0)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/n_cells_region_oi)

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)
    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "real_sampling", real_artificial = "real", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE)

  }

  if (dataset_type == "real_sampling_top1") { # put every cell type just with one region, then for every region, we will filter out the filtered cells that belong to this region  - but keep the real proportions
    metadata = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters)

    filtered_region_celltype = metadata %>% group_by(seurat_clusters, prior_regions) %>% dplyr::count() %>%
      ungroup()

    filtered_region_celltype = filtered_region_celltype %>% group_by(seurat_clusters) %>% top_n(1, n) %>% ungroup() %>% select(-n)
    regions = filtered_region_celltype %>% pull(prior_regions) %>% unique()
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% filter(prior_regions == region_oi) %>% pull(seurat_clusters) %>% unique()
      metadata_region_oi = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters) %>% filter(prior_regions == region_oi & seurat_clusters %in% celltypes_oi)
      n_cells_region_oi = metadata_region_oi %>% pull(cell_id) %>% unique() %>% length()

      present_celltypes_numbers = metadata_region_oi %>% group_by(seurat_clusters) %>% dplyr::count()
      absent_celltypes_numbers = tibble(seurat_clusters = setdiff(celltypes, present_celltypes_numbers$seurat_clusters %>% unique), n = 0)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/n_cells_region_oi)

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)
    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "real_sampling_top1", real_artificial = "real", uniform_diverse = "diverse", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE)

  }


  if(dataset_type == "real_sampling_top1_uniform"){ #  # put every cell type just with one region, then for every region, we will filter out the filtered cells that belong to this region  - make the proportions uniform
    metadata = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters)

    filtered_region_celltype = metadata %>% group_by(seurat_clusters, prior_regions) %>% dplyr::count() %>%
      ungroup()

    filtered_region_celltype = filtered_region_celltype %>% group_by(seurat_clusters) %>% top_n(1, n) %>% ungroup() %>% select(-n)
    regions = filtered_region_celltype %>% pull(prior_regions) %>% unique()
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% filter(prior_regions == region_oi) %>% pull(seurat_clusters) %>% unique()
      other_celltypes = setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble(seurat_clusters = celltypes_oi, n = 1)
      absent_celltypes_numbers = tibble(seurat_clusters = other_celltypes, n = 0)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/length(celltypes_oi))

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)
    })

    names(celltype_frequencies_list) = regions
    properties_df = tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "real_sampling_top1_uniform", real_artificial = "real", uniform_diverse = "uniform", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE)

  }
  if(dataset_type == "real_sampling_top2_overlap"){ # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types is according to real data
    metadata = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters)

    filtered_region_celltype = metadata %>% group_by(seurat_clusters, prior_regions) %>% dplyr::count() %>%
      ungroup()

    filtered_region_celltype = filtered_region_celltype %>% group_by(seurat_clusters) %>% top_n(2, n) %>% ungroup() %>% select(-n)
    regions = filtered_region_celltype %>% pull(prior_regions) %>% unique()
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% filter(prior_regions == region_oi) %>% pull(seurat_clusters) %>% unique()
      metadata_region_oi = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters) %>% filter(prior_regions == region_oi & seurat_clusters %in% celltypes_oi)
      n_cells_region_oi = metadata_region_oi %>% pull(cell_id) %>% unique() %>% length()

      present_celltypes_numbers = metadata_region_oi %>% group_by(seurat_clusters) %>% dplyr::count()
      absent_celltypes_numbers = tibble(seurat_clusters = setdiff(celltypes, present_celltypes_numbers$seurat_clusters %>% unique), n = 0)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/n_cells_region_oi)

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)
    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "real_sampling_top2_overlap", real_artificial = "real", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE)

  }

  if(dataset_type == "real_sampling_top2_overlap_uniform"){ #  # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types is uniform
    metadata = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters)

    filtered_region_celltype = metadata %>% group_by(seurat_clusters, prior_regions) %>% dplyr::count() %>%
      ungroup()

    filtered_region_celltype = filtered_region_celltype %>% group_by(seurat_clusters) %>% top_n(2, n) %>% ungroup() %>% select(-n)
    regions = filtered_region_celltype %>% pull(prior_regions) %>% unique()
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% filter(prior_regions == region_oi) %>% pull(seurat_clusters) %>% unique()
      other_celltypes = setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble(seurat_clusters = celltypes_oi, n = 1)
      absent_celltypes_numbers = tibble(seurat_clusters = other_celltypes, n = 0)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/length(celltypes_oi))

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)
    })

    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "real_sampling_top2_overlap_uniform", real_artificial = "real", uniform_diverse = "uniform", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE)

  }

  if (dataset_type == "real_missing_celltypes_visium"){ # take 4 missing cell types from scRNAseq data that are not used in Visium
    metadata = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters)

    filtered_region_celltype = metadata %>% group_by(seurat_clusters, prior_regions) %>% dplyr::count() %>%
      ungroup()

    filtered_region_celltype = filtered_region_celltype %>% group_by(seurat_clusters) %>% top_n(2, n) %>% ungroup() %>% select(-n)
    regions = filtered_region_celltype %>% pull(prior_regions) %>% unique()
    set.seed(1919)
    leave_out_celltypes = celltypes[sample(seq(length(celltypes)),4)]
    print("cell types that will be left out: ")
    print(leave_out_celltypes)
    celltype_frequencies_list = regions %>% lapply(function(region_oi){
      celltypes_oi = filtered_region_celltype %>% filter(prior_regions == region_oi) %>% pull(seurat_clusters) %>% unique() %>% setdiff(leave_out_celltypes)
      metadata_region_oi = seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, prior_regions, seurat_clusters) %>% filter(prior_regions == region_oi & seurat_clusters %in% celltypes_oi)
      n_cells_region_oi = metadata_region_oi %>% pull(cell_id) %>% unique() %>% length()

      present_celltypes_numbers = metadata_region_oi %>% group_by(seurat_clusters) %>% dplyr::count()
      absent_celltypes_numbers = tibble(seurat_clusters = setdiff(celltypes, present_celltypes_numbers$seurat_clusters %>% unique), n = 0)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/n_cells_region_oi)

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)
      return(celltype_frequencies)
    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "real_missing_celltypes_visium", real_artificial = "real", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = TRUE)

  }

  if(dataset_type == "artificial_sampling_uniform_distinct"){ #   every region consists of different cell types, frequency uniform

    celltype_region_splits =  split( celltypes, sample(nregions, length(celltypes), repl = TRUE) )
    celltype_frequencies_list = seq(nregions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]

      other_celltypes = setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble(seurat_clusters = celltypes_oi, n = 1)
      absent_celltypes_numbers = tibble(seurat_clusters = other_celltypes, n = 0)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/length(celltypes_oi))

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "artificial_sampling_uniform_distinct", real_artificial = "artificial", uniform_diverse = "uniform", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE)

  }
  if(dataset_type == "artificial_sampling_diverse_distinct"){ #  every region consists of different cell types, frequency random
    celltype_region_splits =  split( celltypes, sample(nregions, length(celltypes), repl = TRUE) )
    celltype_frequencies_list = seq(nregions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]

      other_celltypes = setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble(seurat_clusters = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble(seurat_clusters = other_celltypes, n = 0)

      sum_present = sum(present_celltypes_numbers$n)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "artificial_sampling_diverse_distinct", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "distinct", dominant_celltype = FALSE, missing_celltype = FALSE)

  }
  if(dataset_type == "artificial_sampling_uniform_overlap"){ # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types will be uniform

    overlapping_celltypes = celltypes[sample(seq(length(celltypes)),4)] %>% as.character() # and out of these 4, 1 will be chosen per region
    celltype_region_splits =  split( celltypes, sample(nregions, length(celltypes), repl = TRUE) )
    celltype_frequencies_list = seq(nregions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]
      celltypes_oi = union(celltypes_oi, sample(overlapping_celltypes,1))
      other_celltypes = setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble(seurat_clusters = celltypes_oi, n = 1)
      absent_celltypes_numbers = tibble(seurat_clusters = other_celltypes, n = 0)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/length(celltypes_oi))

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "artificial_sampling_uniform_overlap", real_artificial = "artificial", uniform_diverse = "uniform", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE)

  }
  if(dataset_type == "artificial_sampling_diverse_overlap"){ # every region consists of different cell types, but some overlap and replacement is possible - frequency of cell types will be random

    overlapping_celltypes = celltypes[sample(seq(length(celltypes)),4)] %>% as.character() # and out of these 4, 1 will be chosen per region
    celltype_region_splits =  split( celltypes, sample(nregions, length(celltypes), repl = TRUE) )
    celltype_frequencies_list = seq(nregions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]]
      celltypes_oi = union(celltypes_oi, sample(overlapping_celltypes,1))

      other_celltypes = setdiff(celltypes, celltypes_oi)

      present_celltypes_numbers = tibble(seurat_clusters = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble(seurat_clusters = other_celltypes, n = 0)

      sum_present = sum(present_celltypes_numbers$n)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)


    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "artificial_sampling_diverse_overlap", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = FALSE, missing_celltype = FALSE)

  }
  if (dataset_type == "artificial_dominant_celltype_diverse"){ # one cell type has significant higher frequency over all regions - difference in regions determined by other cell types. Some cell types in more than one region.

    dominant_celltype = sample(celltypes,1)
    print("dominant cell type is")
    print(dominant_celltype)

    celltype_region_splits =  split( celltypes, sample(nregions, length(celltypes), repl = TRUE) )
    celltype_frequencies_list = seq(nregions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]] %>% setdiff(dominant_celltype)
      other_celltypes = setdiff(celltypes, celltypes_oi) %>% setdiff(dominant_celltype)

      present_celltypes_numbers = tibble(seurat_clusters = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble(seurat_clusters = other_celltypes, n = 0)
      dominant_celltypes_numbers = tibble(seurat_clusters = dominant_celltype, n = sample(seq(75,100), 1))

      sum_present = sum(present_celltypes_numbers$n) %>% sum(dominant_celltypes_numbers$n)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers) %>% bind_rows(dominant_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)

    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "artificial_dominant_celltype_diverse", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = TRUE, missing_celltype = FALSE)

  }
  if (dataset_type == "artificial_partially_dominant_celltype_diverse"){ # one cell type has significant higher frequency over all regions - difference in regions determined by other cell types. Some cell types in more than one region.

    dominant_celltype = sample(celltypes,1)
    print("dominant cell type is")
    print(dominant_celltype)

    celltype_region_splits =  split( celltypes, sample(nregions, length(celltypes), repl = TRUE) )
    celltype_frequencies_list = seq(nregions) %>% lapply(function(i){
      region_oi = regions[i]
      celltypes_oi = celltype_region_splits[[i]] %>% setdiff(dominant_celltype)
      other_celltypes = setdiff(celltypes, celltypes_oi) %>% setdiff(dominant_celltype)

      present_celltypes_numbers = tibble(seurat_clusters = celltypes_oi, n = sample(10,length(celltypes_oi)))
      absent_celltypes_numbers = tibble(seurat_clusters = other_celltypes, n = 0)
      if(i == length(regions)){
        dominant_celltypes_numbers = tibble(seurat_clusters = dominant_celltype, n = 0)
      } else if (i == length(regions) -1){
        dominant_celltypes_numbers = tibble(seurat_clusters = dominant_celltype, n = sample(seq(1,10), 1))
      } else {
        dominant_celltypes_numbers = tibble(seurat_clusters = dominant_celltype, n = sample(seq(75,100), 1))
      }

      sum_present = sum(present_celltypes_numbers$n) %>% sum(dominant_celltypes_numbers$n)

      celltype_numbers = bind_rows(present_celltypes_numbers, absent_celltypes_numbers) %>% bind_rows(dominant_celltypes_numbers)
      celltype_numbers = celltype_numbers %>% mutate(freq = n/sum_present)

      celltype_frequencies = celltype_numbers %>% pull(freq)
      names(celltype_frequencies) = celltype_numbers %>% pull(seurat_clusters)

      return(celltype_frequencies)

    })
    names(celltype_frequencies_list) = regions
    properties_df =  tibble(dataset_name = paste0(dataset_type, dataset_name),dataset_type = "artificial_partially_dominant_celltype_diverse", real_artificial = "artificial", uniform_diverse = "diverse", distinct_overlap = "overlap", dominant_celltype = TRUE, missing_celltype = FALSE)

  }
  # for in next function, we will need to have the variable in which region is stored if generating 'real synthetic data'
  properties_df = properties_df %>% inner_join(tibble(real_artificial = c("real","artificial"), real_region_var = c(region_var, NA)))

  # Define the number of spots that needs to be generated per region (random number between 100 and 700)
  n_spots_list = regions %>% lapply(function(region_oi){
    return(sample(seq(100,700),1))
  })
  names(n_spots_list) = regions

  celltype_frequencies_list = celltype_frequencies_list[names(n_spots_list)]
  region_assignments = purrr::map2(n_spots_list, celltype_frequencies_list, function(x,y){return(list(n_spots = x, celltype_freq = y))})
  gs = purrr::map2(names(celltype_frequencies_list),celltype_frequencies_list, function(region, celltype_freq){
    tibble(prior_region = region, celltype = names(celltype_freq), freq = celltype_freq, present = celltype_freq > 0)
  }) %>% bind_rows()

  return(list(region_assignments = region_assignments, dataset_properties = properties_df, gold_standard_priorregion = gs))

}
#' @title Get a Seurat object only containing cells of a region (or regions) of interest.
#'
#' @description \code{get_seuratObj_region_seurat} Get a Seurat object only containing cells of a region (or regions) of interest. This filtering is based on spatial annotation of cells as done by `predict_regions_of_cells`, and via the Seurat Label Transfer method.
#' @usage
#' get_seuratObj_region_seurat(region_oi, seurat_obj_scrnaseq, spatial_annotation_assay = "regions", label_transfer_cutoff = 0.25)
#'
#' @inheritParams predict_regions_of_cells
#' @inheritParams get_ligands_receptors_of_region
#' @param spatial_annotation_assay Name of the assay of the Seurat object where the spatial region annotation predictions of cells were stored. Default: 'regions', as resulting from running `predict_regions_of_cells`
#' @param label_transfer_cutoff Cells with a spatial annotation score larger or equal than this cutoff will be kept, and thus considered as belonging to the region of interest. Default = 0.25.
#'
#' @return The Seurat object with the scRNAseq data (seurat_obj_scrnaseq) of cells belonging to the region(s) of interest.
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(Seurat)
#' seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
#' seurat_obj_scrnaseq_central = get_seuratObj_region_seurat(region_oi = "central", seurat_obj_scrnaseq, spatial_annotation_assay = "regions", label_transfer_cutoff = 0.25)
#' }
#'
#' @export
#'
generate_syn_data = function(region_assignment_list, seurat_obj, clust_var, visium_mean = 20000, visium_sd = 7000){

  library(DropletUtils) # For the downsampling
  library(dtplyr) # To use dplyr commands with DT speed

  region_assignments = region_assignment_list[["region_assignments"]]
  dataset_properties = region_assignment_list[["dataset_properties"]]
  real_dataset = dataset_properties$real_artificial == "real"
  gold_standard_priorregion = region_assignment_list[["gold_standard_priorregion"]]

  seurat_obj$seurat_clusters = droplevels(seurat_obj@meta.data[,clust_var] %>% as.factor())
  if(real_dataset == TRUE){
    region_var = dataset_properties$real_region_var
    seurat_obj$prior_regions = droplevels(seurat_obj@meta.data[,region_var] %>% as.factor())
  }

  regions = names(region_assignments)
  synthetic_regions = regions %>% lapply(generate_spots, region_assignments, seurat_obj, visium_mean, visium_sd, real_dataset)

  syn_count_data = synthetic_regions %>% purrr::map("counts") %>% do.call("cbind",.)
  syn_spot_composition = synthetic_regions %>% purrr::map("spot_composition") %>% bind_rows()
  syn_spot_composition_relative = syn_spot_composition %>% .[,-which(colnames(.) %in% c("name","region"))] %>% apply(1, function(x){x/sum(x)}) %>% t() %>% cbind(.,syn_spot_composition[,c("name","region")])
  return(list(counts = syn_count_data, spot_composition = syn_spot_composition, relative_spot_composition = syn_spot_composition_relative, gold_standard_priorregion = gold_standard_priorregion, dataset_properties = dataset_properties))

}
#' @title Get a Seurat object only containing cells of a region (or regions) of interest.
#'
#' @description \code{get_seuratObj_region_seurat} Get a Seurat object only containing cells of a region (or regions) of interest. This filtering is based on spatial annotation of cells as done by `predict_regions_of_cells`, and via the Seurat Label Transfer method.
#' @usage
#' get_seuratObj_region_seurat(region_oi, seurat_obj_scrnaseq, spatial_annotation_assay = "regions", label_transfer_cutoff = 0.25)
#'
#' @inheritParams predict_regions_of_cells
#' @inheritParams get_ligands_receptors_of_region
#' @param spatial_annotation_assay Name of the assay of the Seurat object where the spatial region annotation predictions of cells were stored. Default: 'regions', as resulting from running `predict_regions_of_cells`
#' @param label_transfer_cutoff Cells with a spatial annotation score larger or equal than this cutoff will be kept, and thus considered as belonging to the region of interest. Default = 0.25.
#'
#' @return The Seurat object with the scRNAseq data (seurat_obj_scrnaseq) of cells belonging to the region(s) of interest.
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(Seurat)
#' seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
#' seurat_obj_scrnaseq_central = get_seuratObj_region_seurat(region_oi = "central", seurat_obj_scrnaseq, spatial_annotation_assay = "regions", label_transfer_cutoff = 0.25)
#' }
#'
#' @export
#'
generate_spots = function(region_oi, region_assignments, seurat_obj, visium_mean = 20000, visium_sd = 7000, real_dataset = FALSE){

  count_matrix = as.matrix(seurat_obj[["RNA"]]@counts)
  metadata = seurat_obj@meta.data
  metadata = metadata %>% rownames_to_column("cell_id") %>% as_tibble()

  info_oi = region_assignments[[region_oi]]
  n = info_oi$n_spots
  celltype_freq = info_oi$celltype_freq

  # print(region_oi)
  # print(celltype_freq)
  ds_spots = lapply(1:n, function(i){

    # Select between 2 and 10 cells randomly from the count matrix
    # add here something that uses the cell type frequencies!
    n_cells = sample(x = 2:10, size = 1)

    celltypes_pool = sample(x = names(celltype_freq), size =  n_cells, prob = celltype_freq, replace = T)
    celltypes_pool_table = celltypes_pool %>% table()

    cell_pool = names(celltypes_pool_table) %>% lapply(function(celltype_oi, celltypes_pool_table, metadata, real_dataset){
      n_cells_oi = celltypes_pool_table[celltype_oi]
      # print(n_cells_oi)
      if(real_dataset == FALSE){
        source_cells = metadata %>% filter(seurat_clusters == celltype_oi) %>% pull(cell_id) #
      } else {
        source_cells = metadata %>% filter(seurat_clusters == celltype_oi & prior_regions == region_oi) %>% pull(cell_id) # if dataset is real: pick only cells from possible region! (but then it should be present right...)
      }
      # print(length(source_cells))
      picked_cells = sample(source_cells, size = n_cells_oi, replace = T)
      return(picked_cells)
    }, celltypes_pool_table, metadata, real_dataset) %>% unlist() %>% unique()

    # Create a name for the spot with the necessary info to deconvolute it
    pos = which(colnames(count_matrix) %in% cell_pool)
    tmp_ds = seurat_obj@meta.data[pos,] %>% mutate(weight = 1)
    name_simp = paste(region_oi,'_spot_',i,sep='')

    # table with information of cell types per spots
    spot_ds = suppressMessages(tmp_ds %>%
                                 select(seurat_clusters, weight) %>%
                                 mutate(seurat_clusters = as.character(seurat_clusters)) %>%
                                 # mutate(seurat_clusters = paste('clust_',seurat_clusters,sep = '')) %>%
                                 group_by(seurat_clusters) %>%
                                 summarise(sum_weights = sum(weight)) %>% #
                                 pivot_wider(names_from = seurat_clusters, values_from = sum_weights) %>%
                                 mutate(name = name_simp, region = region_oi))

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

      syn_spot_sparse = downsampleMatrix(Matrix::Matrix(syn_spot,sparse = T), prop = prop_oi)
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
    bind_rows() %>%
    lazy_dt(.) %>%  # Convert to lazy data so it can be avaluated as a data.table
    # mutate_all(~replace_na(.,0)) %>%
    data.frame()

  ds_spots_metadata[is.na(ds_spots_metadata)] = 0
  # change column order so that its progressive
  lev_mod = gsub("[\\+|\\ ]", ".", levels(seurat_obj$seurat_clusters))
  all_cn = c(lev_mod,'name','region')
  if( sum(all_cn %in% colnames(ds_spots_metadata)) == (nlevels(seurat_obj$seurat_clusters)+2) ){
    ds_spots_metadata = ds_spots_metadata[,all_cn]
  } else {
    missing_cols = all_cn[ which( ! all_cn %in% colnames(ds_spots_metadata) ) ]
    ds_spots_metadata[missing_cols] = 0
    ds_spots_metadata = ds_spots_metadata[,all_cn]
  }
  return(list(counts = ds_syn_spots, spot_composition = ds_spots_metadata))

}
