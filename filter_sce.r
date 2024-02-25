filter_sce <- function(){
	
	#get clusters
	sce <- FindClusters(sce, resolution = 2)

	#remove allen cells from object
	sce <- subset(sce, subset = technology == 'tenx')
	total_10x_prefilter <- ncol(sce)
	table_total <-  as.data.frame(table(sce[[]]$categories))

	#remove cluters co-expressing gaba and vglut1
	cluster_exclusion_list <- c('3','7')
	sce <- subset(sce, subset = seurat_clusters %in% cluster_exclusion_list, invert = TRUE)
	cluster_filter_table <- as.data.frame(table(sce[[]]$categories))

	#remove cells with distant neighbors
	sce <- subset(sce, subset = distance < variables$nearest_neighbors_cutoff)
	table_post_neighbor_filter <- as.data.frame(table(sce[[]]$categories))

	#excluded groups with less than 500 cells
	exclusion_list <- c('Astro','Meis2','Peri','Serpinf1','SMC','Sncg','VLMC')
	sce <- subset(sce, subset = categories %in% exclusion_list, invert = TRUE)
	table_after_group_filter <- as.data.frame(table(sce[[]]$categories))

	# format the counts table
	rownames(table_total) <- table_total$Var1
	table_total$Var1 <- NULL
	colnames(table_total) <- c('total')
	table_total$Var1 <- NULL

	rownames(cluster_filter_table) <- cluster_filter_table$Var1
	cluster_filter_table$Var1 <- NULL
	colnames(cluster_filter_table) <- c('post cluster filter')

	rownames(table_post_neighbor_filter) <- table_post_neighbor_filter$Var1
	table_post_neighbor_filter$Var1 <- NULL
	colnames(table_post_neighbor_filter) <- c('post neighbor filter')

	rownames(table_after_group_filter) <- table_after_group_filter$Var
	table_after_group_filter$Var1 <- NULL
	colnames(table_after_group_filter) <- c('post group filter')

	df1 <- transform(merge(table_total,cluster_filter_table,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
	df2 <- transform(merge(df1,table_post_neighbor_filter,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
	df_counts <- transform(merge(df2,table_after_group_filter,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
	df_counts[is.na(df_counts)] <- 0

	# add totals
	df_counts[15,] = c(colSums(df_counts[,1:4]))
	rownames(df_counts)[rownames(df_counts) == "15"] <- "total"

	# write it to a summary file
	stargazer(df_counts, header=FALSE, type='text', summary=FALSE, title="Count summary",digits=1, out = 'count_summary.txt')
	
	# save the processed file
	saveRDS(sce,get_path(variables$stored_variables,'sce.filtered.rds'))
}
