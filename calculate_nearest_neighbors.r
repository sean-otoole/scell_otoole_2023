# Function to calculate nearest neighbors based on iNMF embeddings
calculate_nearest_neighbors <- function() {

  # Subset the SingleCellExperiment (SCE) object into two groups based on technology type
  just.allen <- subset(sce, technology == 'smartseq')  # Cells from the 'smartseq' technology
  our.data <- subset(sce, technology == 'smartseq', invert = TRUE)  # All other data (non-smartseq)

  # Extract iNMF embeddings for both subsets (just.allen and our.data)
  allen.vectors <- Embeddings(object = just.allen[['iNMF']])[,(1:variables$k_value)[variables$removed_dims]]
  our.vectors <- Embeddings(object = our.data[['iNMF']])[,(1:variables$k_value)[variables$removed_dims]]

  # Initialize empty vectors to store results for each iteration
  neighbor.classes <- c()
  euclidean.distances <- c()
  neighbor.categories <- c()
  neighbor.groups <- c()

  # Loop over each row (cell) in the 'our.data' subset to find nearest neighbors
  for (i in 1:dim(our.vectors)[1]) {
    
    # Calculate the Euclidean distance between the current cell's iNMF vector and all 'allen.vectors'
    nearest.neighbor <- sort(sqrt(rowSums((sweep(allen.vectors, 2, our.vectors[i,], '-')^2))))[1:variables$neighbor_value]
    euclidean.distance <- nearest.neighbor[[1]]  # Store the minimum Euclidean distance
    closest.allen.cells <- names(nearest.neighbor)  # Get the names of the nearest cells from 'just.allen'
    
    # Extract the neighborhood class for each of the nearest cells from the Allen dataset
    neighborhood_classes <- just.allen[[]][closest.allen.cells, 4]
    
    # Calculate weights inversely proportional to the Euclidean distance (closer neighbors have higher weights)
    weights <- 1 / nearest.neighbor
    weights.df <- data.frame(neighborhood_classes, weights)
    
    # Sum the weights for each unique neighborhood class
    summed_neighbor_weights <- tapply(weights.df$weights, weights.df$neighborhood_classes, FUN = sum)
    
    # Determine the neighborhood class with the highest summed weight
    neighbor.class <- names(which.max(summed_neighbor_weights))
    class.index <- which(just.allen[[]]$allenCluster == neighbor.class)[1]  # Find the index of the class
    
    # Store the corresponding subclass and group information from the Allen dataset
    neighbor.categories <- rbind(neighbor.categories, just.allen[[]]$allenSubClass[class.index])
    neighbor.groups <- rbind(neighbor.groups, just.allen[[]]$allenClass[class.index])
    neighbor.classes <- rbind(neighbor.classes, neighbor.class)
    euclidean.distances <- rbind(euclidean.distances, euclidean.distance)  # Store the Euclidean distance
  }

  # Ensure the length of 'euclidean.distances' matches the number of columns in 'just.allen'
  for (i in 1:dim(just.allen)[2]) {
    euclidean.distances <- rbind(0, euclidean.distances)  # Padding with zero for alignment
  }

  # Prepend the original Allen dataset information to the results
  neighbor.classes <- c(c(just.allen[[]][, 4][1:length(just.allen[[]][, 4])]), c(neighbor.classes))
  neighbor.categories <- c(c(just.allen[[]][, 5][1:length(just.allen[[]][, 5])]), c(neighbor.categories))
  neighbor.groups <- c(c(just.allen[[]][, 6][1:length(just.allen[[]][, 6])]), c(neighbor.groups))

  # Store the calculated neighbors and distances as new metadata in the SCE object
  sce$neighbors <- c(neighbor.classes)  # Nearest neighbor class for each cell
  sce$distance <- c(euclidean.distances)  # Euclidean distance to the nearest neighbor
  sce$categories <- c(neighbor.categories)  # Subcategory of the nearest neighbor
  sce$groups <- c(neighbor.groups)  # Group of the nearest neighbor

  # Save the updated SCE object with the calculated nearest neighbor information
  saveRDS(sce, get_path(variables$stored_variables, 'sce.post.nearest.neighbors.rds'))
}
