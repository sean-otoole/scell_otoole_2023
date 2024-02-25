get_variables <- function() {
	variable_list <- list("mito_max" = 20, "mito_min" = 0, "FDR_thresh" = 0.05, "cosine_threshold" = 0,
	"initial_features" = 2000, "exclude_features" = TRUE, "positive_markers" = TRUE, 
	"k_value" = 12, "lambda_value" = 10,"neighbor_values" = 10,
	"removed_dims" = c(-16), "nearest_neighbors_cutoff" = 0.002,
	"data_repository" = 'data_repository/',
	"stored_variables" = 'stored_variables/',
	"figures" = 'code_generated_figures/',
	"FixAxes" = 'FALSE',
	"BootStrapIters" = 100000)
	return(variable_list) 
}

get_variables()
