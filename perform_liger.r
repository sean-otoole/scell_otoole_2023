perform_liger <- function() {
sce <- NormalizeData(sce)
VariableFeatures(sce) <- features
sce <- ScaleData(sce, split.by = "facs_sample", do.center = FALSE, features = features)
sce <- RunOptimizeALS(sce, k = variables$k_value, lambda = variables$lambda_value, split.by = "facs_sample",ref_dataset = 'allen')
sce <- RunQuantileNorm(sce, split.by = "facs_sample",ref_dataset = 'allen')
sce <- FindNeighbors(sce, reduction = "iNMF", dims = ((1:ncol(sce[["iNMF"]]))[variables$removed_dims]))
sce <- FindClusters(sce, resolution = 0.3)
sce <- RunUMAP(sce, dims = (1:ncol(sce[["iNMF"]])), reduction = "iNMF", n.neighbors = 30, min.dist = 0.4, spread = 0.5,
              repulsion.strength = 0.2)
saveRDS(sce,get_path(variables$stored_variables,'sce.post.liger.rds'))
}
