calculate_nearest_neighbors <- function(){

	just.allen <- subset(sce, technology == 'smartseq')
	our.data <- subset(sce, technology == 'smartseq',invert = TRUE)

	allen.vectors <- Embeddings(object = just.allen[['iNMF']])[,(1:variables$k_value)[variables$removed_dims]]
	our.vectors <- Embeddings(object = our.data[['iNMF']])[,(1:variables$k_value)[variables$removed_dims]]

	neighbor.classes <- c()
	euclidean.distances <- c()
	neighbor.categories <- c()
	neighbor.groups <- c()

	for (i in 1:dim(our.vectors)[1]){
    		nearest.neighbor <- sort(sqrt(rowSums((sweep(allen.vectors,2,our.vectors[i,],'-')^2))))[1:variables$neighbor_value]
    		euclidean.distance <- nearest.neighbor[[1]]
    		closest.allen.cells <- names(nearest.neighbor)
    		neighborhood_classes <- just.allen[[]][closest.allen.cells,4]
    		weights <- 1/nearest.neighbor
    		weights.df <- data.frame(neighborhood_classes, weights)
    		summed_neighbor_weights <- tapply(weights.df$weights, weights.df$neighborhood_classes, FUN=sum)
    		neighbor.class <- names(which.max(summed_neighbor_weights))
    		class.index <- which(just.allen[[]]$allenCluster == neighbor.class)[1]
    		neighbor.categories <- rbind(neighbor.categories,just.allen[[]]$allenSubClass[class.index])
    		neighbor.groups <- rbind(neighbor.groups,just.allen[[]]$allenClass[class.index])
    		neighbor.classes <- rbind(neighbor.classes,neighbor.class)
    		euclidean.distances <- rbind(euclidean.distances,euclidean.distance)
	}

	for (i in 1:dim(just.allen)[2]){
    		euclidean.distances <- rbind(0, euclidean.distances)
	}

	neighbor.classes <- c(c(just.allen[[]][,4][1:length(just.allen[[]][,4])]), c(neighbor.classes))
	neighbor.categories <- c(c(just.allen[[]][,5][1:length(just.allen[[]][,5])]), c(neighbor.categories))
	neighbor.groups <- c(c(just.allen[[]][,6][1:length(just.allen[[]][,6])]), c(neighbor.groups))

	sce$neighbors <- c(neighbor.classes)
	sce$distance <- c(euclidean.distances)
	sce$categories <- c(neighbor.categories)
	sce$groups <- c(neighbor.groups)

	saveRDS(sce,get_path(variables$stored_variables,'sce.post.nearest.neighbors.rds'))
}
