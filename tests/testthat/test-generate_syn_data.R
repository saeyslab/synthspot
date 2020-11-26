context("generating synthetic data works for each 'tissue type'")
test_that("generating synthetic data works for each 'tissue type' of 'real' class", {


  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real", clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real_top1", clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real_top1_uniform", clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real_top2_overlap", clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real_top2_overlap_uniform", clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")


})


test_that("generating synthetic data works for 'tissue types' with missing cell types", {

  set.seed(1919) # there might be a problem here with the toy data and missing cell types...
  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real_missing_celltypes_visium", clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  set.seed(1919) # there might be a problem here with the toy data and missing cell types...
  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_missing_celltypes_visium", clust_var = "subclass", region_var = NULL , n_regions = 11,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

})

test_that("generating synthetic data works for each 'tissue type' of 'artificial' class", {


  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_uniform_distinct", clust_var = "subclass", region_var = NULL , n_regions = 5,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_diverse_distinct", clust_var = "subclass", region_var = NULL , n_regions = 5,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_uniform_overlap", clust_var = "subclass", region_var = NULL , n_regions = 5,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_diverse_overlap", clust_var = "subclass", region_var = NULL , n_regions = 5,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_dominant_celltype_diverse", clust_var = "subclass", region_var = NULL , n_regions = 5,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_partially_dominant_celltype_diverse", clust_var = "subclass", region_var = NULL , n_regions = 5,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_uniform_distinct", clust_var = "subclass", region_var = NULL , n_regions = 3,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000)
  expect_type(object = synthetic_visium_data, type = "list")

  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_uniform_distinct", clust_var = "subclass", region_var = NULL , n_regions = 3,
                                                    n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000, add_mock_region = TRUE)
  expect_type(object = synthetic_visium_data, type = "list")




})
test_that("input checks work properly", {

  expect_error(synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_uniform_distinct", clust_var = "celltype_class", region_var = NULL , n_regions = 5,
                                                                 n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000))
  expect_error(synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_uniform_distinct", clust_var = 1, region_var = NULL , n_regions = 5,
                                                                 n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000))
  expect_error(synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real", clust_var = "subclass", region_var = NULL , n_regions = 1,
                                                                 n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000))
  expect_error(synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "real", clust_var = "subclass", region_var = 1 , n_regions = NULL,
                                                                 n_spots_min = 5, n_spots_max = 20, visium_mean = 20000, visium_sd = 5000))

})
