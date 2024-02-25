# preprocessing

preprocess <- function(input_directory){

    #Read in the data

    sce <- read10xCounts(input_directory, col.names = TRUE, type = c("auto"), version = ('auto'))

    #set proper gene names

    rowDataInfo <- (rowData(sce))
    Symbols <- (rowDataInfo[,2])
    duplicateNumber <- 1
    while (length(Symbols[duplicated(Symbols)])> 0)
    {
        duplicateString <- toString(duplicateNumber)
        duplicate_symbols <- Symbols[duplicated(Symbols)]
        replacements <- paste0(Symbols[duplicated(Symbols)],'-',duplicateString)
        redundant_indexes <- match(duplicate_symbols,Symbols)
        Symbols[redundant_indexes] <- replacements
        row.names(sce) <- Symbols
        duplicateNumber <- duplicateNumber + 1
    }
    
    #droplet calling and quality filtering

    set.seed(100)
    sum_list <- (unname(colSums(counts(sce))))
    sum_list <- sum_list[sum_list > 10]
    lower_limit <- median(sum_list)
    e.out <- emptyDrops(counts(sce), BPPARAM = MulticoreParam(), lower = lower_limit)
    sce <- sce[,which(e.out$FDR <= variables$FDR_thresh)]
    is.mito <- grep("^mt-", rowData(sce)$Symbol)
    sce <- addPerCellQC(sce, subsets=list(MT=is.mito))
    high.mito <- sce$subsets_MT_percent > variables$mito_max
    min.mito <- sce$subsets_MT_percent < variables$mito_min
    low.lib <- isOutlier(sce$total, log=TRUE, nmads=3, type="lower")
    discard <- low.lib | high.mito | min.mito
    sce <- sce[,!discard]
    lower_bounds <- quantile(colData(sce)$total, probs = 0.3)
    low_read_cells <- sce$total < lower_bounds
    sce <- sce[,!low_read_cells]
    
    #adding in normalized counts and log counts
    sf <- 2^rnorm(ncol(sce))
    sf <- sf/mean(sf)
    normcounts(sce) <- t(t(counts(sce))/sf)
    logcounts(sce) <- log2(normcounts(sce)+1)

    #created comptabile coldata for seurat object, will be used for vars.to.regress in liger later one

    colData(sce)$'percent.mt' <- sce$subsets_MT_percent
    colData(sce)$'nCount_RNA' <- sce$total	

    return(sce)
}

pathConstructor <- function(sample){
    beginning <- 'data_repository/single_cell_data/'
    end <- '/outs/raw_feature_bc_matrix'
    new_path <- paste(beginning, sample, end, sep = '')
    return(new_path)
}

# preprocess the data

preprocess_pipeline <- function(){

sample_paths <- lapply(sample_list, pathConstructor)

sce <- lapply(sample_paths, preprocess)   # converts the sce objects to a list

# add the batch names for liger

for (sample in 1:length(sce)){
    current_sample <- unique(colData(sce[[sample]])$Sample)
    res <- str_match(current_sample, "single_cell_data/\\s*(.*?)\\s*/outs")
    current_sample <- (res[,2])
    batch <- gsub('[0-9]+', '', current_sample)
    colData(sce[[sample]]) <- cbind(colData(sce[[sample]]), Batch = batch)
    colData(sce[[sample]]) <- cbind(colData(sce[[sample]]), facs_sample = current_sample)
}

saveRDS(sce,get_path(variables$stored_variables,'sce.preprocess.rds'))

}