generate_allen_data <- function(){
allen_path = 'data_repository/single_cell_data/allen_data/mouse_VISp_2018-06-14_exon-matrix.csv'
allen.rna <- as.sparse(read.csv(file = allen_path, sep = ",", header = TRUE, row.names = 1))
gene_path <- 'data_repository/single_cell_data/allen_data/mouse_VISp_2018-06-14_genes-rows.csv'
gene.names <- read.csv(file = gene_path, header = TRUE, sep = ',' )
gene_symbol <- gene.names[ ,'gene_symbol']
rownames(allen.rna) <- gene_symbol
meta_path <-  'data_repository/single_cell_data/allen_data/mouse_VISp_2018-06-14_samples-columns.csv'
sample.info <- read.csv(file = meta_path, header = TRUE, sep = ',' )
allen_cluster <- sample.info[ ,'cluster']
allen_cluster <- sapply(allen_cluster, toString)
allen_sub_class <- sample.info[ ,'subclass']
allen_sub_class <- sapply(allen_sub_class, toString)
allen_class <- sample.info[ ,'class']
allen_class <- sapply(allen_class, toString)
allen_region <- sample.info[ ,'brain_region']
allen_region <- sapply(allen_region, toString)
allen_subregion <- sample.info[ ,'brain_subregion']
allen_subregion <- sapply(allen_subregion, toString)
allen.cells <- CreateSeuratObject(counts = allen.rna)
allen.cells@meta.data$allenCluster <- allen_cluster
allen.cells@meta.data$allenSubClass <- allen_sub_class
allen.cells@meta.data$allenClass <- allen_class
allen.cells@meta.data$allenRegion <- allen_region
allen.cells@meta.data$allenSubregion <- allen_subregion
allen.cells@meta.data$sample <- 'allen'
allen.cells@meta.data$facs_sample <- 'allen'
allen.cells@meta.data$Batch <- 'allen'
allen.cells@meta.data$technology <- 'smartseq'
allen.cells@meta.data$percent.mt <- sample.info$percent_mt_exon_reads
allen.cells@meta.data$percent.rb <- sample.info$percent_rrna_reads
allen.cells@meta.data$percent.virus <- 0
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
layers.list <- c("L1","L2/3","L4","L1-L4","L1-L2/3","L2/3-L4")
cluster.list <- c('L2/3 IT VISp Adamts2','L2/3 IT VISp Rrad','L2/3 IT VISp Agmat','L4 IT VISp Rspo1')
exclusion.list <- c('High Intronic','No Class', 'Low Quality','Batch Grouping','Doublet','L5 IT','L6 IT','L5 PT','Endo','CR','NP')
allen.cells <- subset(x = allen.cells, subset = allenCluster %in% cluster.list 
               | allenSubregion %in% layers.list
               & allenSubClass %not in% exclusion.list)
allen.cells <- NormalizeData(allen.cells)
allen.cells <- FindVariableFeatures(allen.cells, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
allen.cells <- ScaleData(allen.cells)
allen.cells <- RunPCA(allen.cells, verbose = FALSE)
allen.cells <- RunUMAP(allen.cells, reduction = "pca", dims = 1:25, verbose = FALSE)
Idents(object = allen.cells) <- 'allenSubClass'
saveRDS(allen.cells, 'stored_variables/allen.cells.rds')
}
