mergeListObject <- function(input_list){
    for (i in 1:length(input_list)){
        if(i == 1){
            combined_sce <- input_list[[i]]
        } else {
            combined_sce <- cbind(combined_sce, input_list[[i]])
        }
    }
    return(combined_sce)
}
                            
splitSamplesToSeurat <- function(sce){
    i <- 1
    sample_list <- unique(sce$Batch)
    object_list <- list()
    for (name in sample_list){
        included <- colData(sce)$Batch == name 
        sce_subset <- sce[,included]
        sce_seurat <- CreateSeuratObject(counts = counts(sce_subset),logcounts = logcounts(sce_subset), meta.data = as.data.frame(colData(sce_subset)))
        sce_seurat@meta.data$sample <- name
        object_list[[i]] <- sce_seurat
        i <- i + 1
    }
    return(object_list)
}


mergeThem <- function(){

sce <- mergeListObject(sce)     # list is converted to one object

metadata(sce) <- list()   # wipe the metadata, interferes with conversion to Seurat
unique_set <- make.unique(colnames(sce))   #make all cell or colnames unique across sets
colnames(sce) <- unique_set
sce <- splitSamplesToSeurat(sce)

sce <- merge(allen.cells, y = sce)

sce <- NormalizeData(sce) # normalization is neccesary for downstream processes

saveRDS(sce,get_path(variables$stored_variables,'sce.post.merge.rds'))


}
