perform_scmap <- function() {

allen <- readRDS(get_path(variables$stored_variables,'allen.cells.rds'))
Idents(object = allen) <- "allenSubClass"
allen.l23 <- subset(x = allen, idents = c('L2/3 IT'))
allen.sce <- as.SingleCellExperiment(allen.l23)
rowData(allen.sce)$feature_symbol <- rownames(allen.sce)
allen.sce <- allen.sce[!duplicated(rownames(allen.sce)), ]

tenx.l23 <- subset(sce, subset = categories == 'L2/3 IT')

allen.sce <- selectFeatures(allen.sce, n_features = 250, suppress_plot = TRUE)
allen.sce <- indexCluster(allen.sce, cluster_col = 'allenCluster')

names_for_sct <- rownames(allen.sce)[rowData(allen.sce)$scmap_features]

tenx.sub <- subset(x = sce, features = names_for_sct)
tenx.sub <- SCTransform(tenx.sub, vars.to.regress = c("percent.mt","chemistry","facs_sample","nFeature_RNA"), verbose = FALSE)
Idents(tenx.sub) <- 'categories'
tenx.l23.sub <- subset(x = tenx.sub, idents = 'L2/3 IT')

tenx.sce <- as.SingleCellExperiment(tenx.l23.sub, assay = 'SCT')
counts(tenx.sce) <- tenx.l23.sub[['SCT']]@counts
logcounts(tenx.sce) <- as(log2(counts(tenx.sce)+1),'sparseMatrix')

rowData(tenx.sce)$feature_symbol <- rownames(tenx.sce)
tenx.sce <- tenx.sce[!duplicated(rownames(tenx.sce)), ]

#heatmap(as.matrix(metadata(allen.sce)$scmap_cluster_index))  

scmapCluster_results <- scmapCluster(
  projection = tenx.sce,threshold = variables$cosine_threshold, 
  index_list = list(
    allen = metadata(allen.sce)$scmap_cluster_index
  )
)

saveRDS(scmapCluster_results,get_path(variables$stored_variables,'scmapCluster_results.rds'))

tenx.l23@meta.data$scmap_cluster_index <- scmapCluster_results$scmap_cluster_labs

saveRDS(tenx.l23,get_path(variables$stored_variables,'l23.10x.post.scmap.rds'))
}
