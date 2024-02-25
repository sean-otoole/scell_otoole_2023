library('ggplot2')
library('scmap')
library('Seurat')
library('lsa')
library('SummarizedExperiment')
library('SingleCellExperiment')
library('proxyC')
library('ggpubr')
library('grDevices')
library('zoo')
library('abind')
library('ggpubr')
library('grid')
library('gplots')

source('get_path.r')
source('single_cell_variables.r')

variables <- get_variables()

allen <- readRDS(get_path(variables$stored_variables,'allen.cells.rds'))
sce <- readRDS(get_path(variables$stored_variables,'sce.filtered.rds'))

# subset the allen class
Idents(object = allen) <- "allenSubClass"
allen.l23 <- subset(x = allen, idents = c('L2/3 IT'))
allen.sce <- as.SingleCellExperiment(allen.l23)
rowData(allen.sce)$feature_symbol <- rownames(allen.sce)
allen.sce <- allen.sce[!duplicated(rownames(allen.sce)), ]

# subset the tenx experiment and select for l2/3
all.tenx <- subset(x = sce, subset = technology == 'tenx')
tenx.l23 <- subset(all.tenx, subset = categories == 'L2/3 IT')

# select features and index clutser names
allen.sce <- selectFeatures(allen.sce, n_features = 250, suppress_plot = TRUE)
allen.sce <- indexCluster(allen.sce, cluster_col = 'allenCluster')

# get the features that have been selected in alen dataset
names_for_sct <- rownames(allen.sce)[rowData(allen.sce)$scmap_features]

# subset the tenx dataset only for features that have been selected for
# adjust the counts of the features based on mitocondira, chemistry facs_sample and nFeature_RNA
# then subset for layer 2/3
tenx.sub <- subset(x = all.tenx, features = names_for_sct)
tenx.sub <- SCTransform(tenx.sub, vars.to.regress = c("percent.mt","chemistry","facs_sample","nFeature_RNA"), verbose = FALSE)
Idents(tenx.sub) <- 'categories'
tenx.l23.sub <- subset(x = tenx.sub, idents = 'L2/3 IT')

# convert to a single cell experiment assign SCT transformed counts to main counts the convert to logcounts (base 2)
tenx.sce <- as.SingleCellExperiment(tenx.l23.sub, assay = 'SCT')
counts(tenx.sce) <- tenx.l23.sub[['SCT']]@counts
logcounts(tenx.sce) <- as(log2(counts(tenx.sce)+1),'sparseMatrix')

# assign row names
rowData(tenx.sce)$feature_symbol <- rownames(tenx.sce)
tenx.sce <- tenx.sce[!duplicated(rownames(tenx.sce)), ]

# perform scmap and assign results
scmapCluster_results <- scmapCluster(
  projection = tenx.sce,threshold = 0, 
  index_list = list(
    allen = metadata(allen.sce)$scmap_cluster_index
  )
)

colData(tenx.sce)$scmap_cluster_index <- scmapCluster_results$scmap_cluster_labs
colData(tenx.sce)$scmap_cluster_score <- scmapCluster_results$scmap_cluster_siml

# corrleation heatmap code, also scaling within gene expression values
tenx.names <- rownames(logcounts(tenx.sce))
allen.exp.vals <- logcounts(allen.sce[tenx.names,])
allen.exp.vals <- t(scale(t(allen.exp.vals)))
tenx.exp.vals <- logcounts(tenx.sce)
tenx.exp.vals <- t(scale(t(tenx.exp.vals)))

# reorder based on identity
tenx_sort_order <- order(colData(tenx.sce)$scmap_cluster_index)
tenx.exp.vals <- tenx.exp.vals[,tenx_sort_order]
allen_sort_order <- order(colData(allen.sce)$allenCluster)
allen.exp.vals <- allen.exp.vals[,allen_sort_order]

num_columns_to_average = 20

# Read in the matrix
matrix <- as.matrix(tenx.exp.vals)

# Calculate the number of sets of columns to average
num_sets = ncol(matrix) / num_columns_to_average

# Initialize a new matrix to hold the averaged columns
averaged_matrix <- matrix(NA, nrow(matrix), num_sets)

# Loop through each set of columns and average them
for(i in 1:num_sets) {
  start_col = (i - 1) * num_columns_to_average + 1
  end_col = start_col + num_columns_to_average - 1
  averaged_matrix[, i] <- rowMeans(matrix[, start_col:end_col])
}

# View the new matrix with the averaged columns
tenx.exp.vals.avg <- averaged_matrix

# get correlation matrix
tenx_v_allen_corr <- cor(as.matrix(allen.exp.vals),tenx.exp.vals.avg,method = 'spearman')

#change orientation of correlations matrix for image function
tenx_v_allen_corr <- apply(tenx_v_allen_corr, 2, rev)
tenx_v_allen_corr <- t(tenx_v_allen_corr)

colors <- colorRampPalette(c("darkgrey","white","darkred"))(256)
rotate <- function(x) t(apply(x, 2, rev))

fig3c_path <- get_path(variables$figures,"figure_3c.pdf")
pdf(fig3c_path )
heatmap.2(x = apply(t(tenx_v_allen_corr),2,rev), Rowv = NA, Colv = NA,scale='none', col = colors,
          labRow=FALSE, labCol = FALSE, breaks = seq(-0.2, 0.2, length.out = 257), dendrogram = 'none',
          trace = 'none', density.info="none",
	  main = sprintf("10x count:%d, allen count:%d", ncol(tenx.exp.vals),ncol(allen.exp.vals)))
dev.off()

# construct matrices for current assigned 10x groups
adam.tenx.matrix <- logcounts(tenx.sce[,colData(tenx.sce)$scmap_cluster_index == 'L2/3 IT VISp Adamts2'])
agmat.tenx.matrix <- logcounts(tenx.sce[,colData(tenx.sce)$scmap_cluster_index == 'L2/3 IT VISp Agmat'])
rrad.tenx.matrix <- logcounts(tenx.sce[,colData(tenx.sce)$scmap_cluster_index == 'L2/3 IT VISp Rrad'])

# construct the average expression vectors for allen types
adam.allen.exp <- rowMeans(logcounts(allen.sce[tenx.names,colData(allen.sce)$allenCluster == 'L2/3 IT VISp Adamts2']))
agmat.allen.exp <- rowMeans(logcounts(allen.sce[tenx.names,colData(allen.sce)$allenCluster == 'L2/3 IT VISp Agmat']))
rrad.allen.exp <- rowMeans(logcounts(allen.sce[tenx.names,colData(allen.sce)$allenCluster == 'L2/3 IT VISp Rrad']))

# get the average correlation for each
adam_corr_assigned <- mean(cor(as.matrix(adam.tenx.matrix), adam.allen.exp,method = 'spearman'))
agmat_corr_assigned <- mean(cor(as.matrix(agmat.tenx.matrix), agmat.allen.exp,method = 'spearman'))
rrad_corr_assigned <- mean(cor(as.matrix(rrad.tenx.matrix), rrad.allen.exp,method = 'spearman'))

# get counts of each type for boostrapping
adam_count <- ncol(adam.tenx.matrix)
agmat_count <- ncol(agmat.tenx.matrix)
rrad_count <- ncol(rrad.tenx.matrix)

# construct the bootstrap matrix
a1 <- cor(as.matrix(adam.tenx.matrix), adam.allen.exp,method = 'spearman')
b1 <- cor(as.matrix(agmat.tenx.matrix), adam.allen.exp,method = 'spearman')
c1 <- cor(as.matrix(rrad.tenx.matrix), adam.allen.exp,method = 'spearman')
col_1 <- abind(a1,b1,c1, along=1)
a2 <- cor(as.matrix(adam.tenx.matrix), agmat.allen.exp,method = 'spearman')
b2 <- cor(as.matrix(agmat.tenx.matrix), agmat.allen.exp,method = 'spearman')
c2 <- cor(as.matrix(rrad.tenx.matrix), agmat.allen.exp,method = 'spearman')
col_2 <- abind(a2,b2,c2, along=1)
a3 <- cor(as.matrix(adam.tenx.matrix), rrad.allen.exp,method = 'spearman')
b3 <- cor(as.matrix(agmat.tenx.matrix), rrad.allen.exp,method = 'spearman')
c3 <- cor(as.matrix(rrad.tenx.matrix), rrad.allen.exp,method = 'spearman')
col_3 <- abind(a3,b3,c3, along=1)
bootstrap_matrix <- cbind(col_1,col_2,col_3)

## randomly samples columns porportionally to cell counts and returns the average correlation value

start_time <- Sys.time()

original_assignment_vector <- c(rep(1,adam_count),rep(2,agmat_count),rep(3,rrad_count))
loop_count <- variables$BootStrapIters
bootstrap_vector_adam <- vector(mode = "double", length = loop_count)
bootstrap_vector_agmat <- vector(mode = "double", length = loop_count)
bootstrap_vector_rrad <- vector(mode = "double", length = loop_count)


for (x in 1:loop_count) {
set.seed(x)
random_assignment_vector <- sample(original_assignment_vector) 
bootstrap_vector_adam[[x]] <- mean(bootstrap_matrix[which(random_assignment_vector == 1),1])
bootstrap_vector_agmat[[x]] <- mean(bootstrap_matrix[which(random_assignment_vector == 2),2])
bootstrap_vector_rrad[[x]] <- mean(bootstrap_matrix[which(random_assignment_vector == 3),3])
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)

return_p <- function(bootsrap_vector,assigned_corr) {
dist_mean <- mean(bootsrap_vector)
current_sd <- sd(bootsrap_vector)
current_mean <- assigned_corr
current_z_score <- (current_mean - dist_mean)/current_sd
print(current_z_score)
}

random_corr_df_adam <- data.frame(random_corr=bootstrap_vector_adam)
random_corr_df_agmat <- data.frame(random_corr=bootstrap_vector_agmat)
random_corr_df_rrad <- data.frame(random_corr=bootstrap_vector_rrad)

adam_plot <- ggplot(random_corr_df_adam,aes(x=random_corr)) + geom_histogram(color="black", fill="gray", binwidth=0.00001) + 
xlim(0.2,0.4) +
theme_classic() + geom_vline(aes(xintercept = adam_corr_assigned), colour="orange") +
ggtitle('bootstrap adamts2')

agmat_plot <- ggplot(random_corr_df_agmat,aes(x=random_corr)) + geom_histogram(color="black", fill="gray", binwidth=0.00001) + 
xlim(0.275,0.325) + theme_classic() + 
geom_vline(aes(xintercept = agmat_corr_assigned), colour="black") +
ggtitle('bootstrap agmat')

rrad_plot <- ggplot(random_corr_df_rrad,aes(x=random_corr)) + geom_histogram(color="black", fill="gray", binwidth=0.00001) + 
xlim(0.2,0.3) + theme_classic() + 
geom_vline(aes(xintercept = rrad_corr_assigned), colour="blue") +
ggtitle('bootstrap rrad')

# textGrob(paste('z_score:',return_p(bootstrap_vector_adam,adam_corr_assigned)))

adam_z <- paste('z_score:',return_p(bootstrap_vector_adam,adam_corr_assigned))
agmat_z <- paste('z_score:',return_p(bootstrap_vector_agmat,agmat_corr_assigned))
rrad_z <- paste('z_score:',return_p(bootstrap_vector_rrad,rrad_corr_assigned))

adam_grob <- grobTree(textGrob(adam_z, x=0.1,  y=0.95, hjust=0,
  gp=gpar(col="red", fontsize=13)))
agmat_grob <- grobTree(textGrob(agmat_z, x=0.1,  y=0.95, hjust=0,
  gp=gpar(col="red", fontsize=13)))
rrad_grob <- grobTree(textGrob(rrad_z, x=0.1,  y=0.95, hjust=0,
  gp=gpar(col="red", fontsize=13)))

adam_bootstrap_path <- get_path(variables$figures,"adam_boostrap.pdf")
ggsave(adam_bootstrap_path,adam_plot)

agmat_bootstrap_path <- get_path(variables$figures,"agmat_boostrap.pdf")
ggsave(agmat_bootstrap_path,agmat_plot)

rrad_bootstrap_path <- get_path(variables$figures,"rrad_boostrap.pdf")
ggsave(rrad_bootstrap_path,rrad_plot)

