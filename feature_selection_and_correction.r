feature_selection_and_correction <- function(){

# split the data
Idents(sce) <- 'facs_sample'
allen <- subset(x = sce, idents = 'allen')fr
tenx <- subset(x = sce, idents = 'allen', invert = TRUE)

# Add chemistry information
tenx.v2 <- subset(x = tenx, subset = Batch %in% c('runa','visa'), invert = FALSE)
tenx.v3 <- subset(x = tenx, subset = Batch %in% c('runa','visa'), invert = TRUE)
tenx.v2$chemistry <- 'v2'
tenx.v3$chemistry <- 'v3'
tenx <- merge(tenx.v2, y = tenx.v3)

#scale the counts and reassign
Idents(tenx) <- 'facs_sample'
tenx <- SCTransform(tenx, vars.to.regress = c("percent.mt","chemistry","facs_sample","nFeature_RNA"), verbose = FALSE, variable.features.n = variables$initial_features)
sctFeatures <- VariableFeatures(tenx)
sctcounts <- tenx[["SCT"]]@counts
sctcounts <- CreateAssayObject(counts = sctcounts)
tenx[["RNA"]] <- sctcounts
DefaultAssay(tenx) <- 'RNA'
tenx$technology <- 'tenx'
all.data <- merge(allen, y = tenx)
saveRDS(all.data,get_path(variables$stored_variables,'preliger.rds'))

#excluding markers associated with cell state (low read count, high mito pct. etc

variables$exclude_features
markers.to.exclude <- readRDS('markers.to.exclude.rds'))
sct.features.filtered <- sctFeatures[!(sctFeatures %in% markers.to.exclude)]
saveRDS(sct.features.filtered,get_path(variables$stored_variables,'filtered.features.rds'))


if (variables$exclude_features == TRUE) {
markers.to.exclude <- readRDS('markers.to.exclude.rds'))
sce <- readRDS(get_path(variables$stored_variables,'preliger.rds'))

sct.features.filtered <- sctFeatures[!(sctFeatures %in% markers.to.exclude)]
saveRDS(sct.features.filtered,get_path(variables$stored_variables,'filtered.features.rds'))
} else {
saveRDS(sctFeatures,get_path(variables$stored_variables,'filtered.features.rds'))
}

}
