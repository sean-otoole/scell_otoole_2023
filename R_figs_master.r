# Suppress package startup messages for a cleaner output
suppressPackageStartupMessages(library('crayon'))

# Inform the user that packages are being loaded
cat(red("\nLoading packages, this will take a few seconds\n"))

# List of required packages for the analysis
Packages <- c("DropletUtils", "scater", "TxDb.Mmusculus.UCSC.mm10.ensGene", 'lsa', 'zoo', 'grid',
              'AnnotationHub', 'scRNAseq', 'scran', 'Matrix', 'BiocParallel', 'scales', 'ggpubr',
              'BiocSingular', 'dplyr', 'edgeR', 'pheatmap', 'extrafont', 'SingleCellExperiment',
              'Seurat', 'stringr', 'SeuratWrappers', 'pdftools', 'tidyverse', 'gplots', 'gridExtra',
              'dendextend', 'RColorBrewer', 'ggplot2', 'png', 'liger', 'rliger', 'scmap', 'tidyr', 'stargazer', 'qdapRegex')

# Load each package silently without startup messages
suppressPackageStartupMessages(invisible(lapply(Packages, library, character.only = TRUE)))

# Inform the user that package loading is complete
cat(red("\nFinished with package loading\n"))

# Source external R scripts (functions for different parts of the analysis pipeline)
source('get_path.r')                          # Function to retrieve file paths
source('single_cell_variables.r')             # Variables related to the single cell analysis
source('generate_allen_data.r')               # Function to generate allen data
source('preprocess.r')                        # Preprocessing steps for data
source('merge_and_split.r')                   # Merge and split data for analysis
source('feature_selection_and_correction.r')  # Feature selection and correction
source('perform_liger.r')                     # Perform LIGER analysis
source('calculate_nearest_neighbors.r')       # Calculate nearest neighbors
source('filter_sce.r')                        # Filter SingleCellExperiment (SCE) object
source('l23_scmap_assignment.r')              # SCmap subtype assignment for Layer 2/3
source('get_the_figs.r')                      # Function to generate the figures

# Load variables needed for the analysis
variables <- get_variables()

# Define the list of samples for the analysis
sample_list <- c('vis1a','vis1b','vis2a','vis2b','vis3a','vis3b','vis4b',
                 'mm1a','mm1b','mm2a','mm2b','mm3a','mm3b','mm4a','mm4b',
                 'run1a','run1b','run2a','run2b','run3a','run3b','run4a','run4b')

# Display the main options to the user
cat("\n")
cat(red("Select a start point (indicate by # value for option): \n
        \t 1.Generate the allen reference
        \t 2.Pre-process the data
        \t 3.Merge and split
        \t 4.GLM correction and feature selection
        \t 5.Run the liger pipeline
        \t 6.Execute the nearest neighbor code and 
        \t 7.filter clusters/distant neighbors
        \t 8.scmap layer2/3 subtype assignment
        \t 9.generate all R related figures"))

# Read user input for the selected analysis step
cat("\n")
user.selection <- readLines(con = "stdin", 1)

# Step 1: Generate the Allen data
if (user.selection == '1') {
  cat(red("\nGenerating the allen data"))
  cat("\n")
  generate_allen_data()  # Generate Allen data
  user.selection <- '2'  # Move to the next step
}

# Step 2: Pre-process the data
if (user.selection == '2') {
  cat(red("\nPre-processing the data"))
  cat("\n")
  preprocess_pipeline()  # Execute the preprocessing pipeline
  user.selection <- '3'  # Move to the next step
} 

# Step 3: Merge and split the data
if (user.selection == '3') {
  cat(red("\nPerforming merge and split\n"))
  cat("\n")
  allen.cells <- readRDS(get_path(variables$stored_variables, 'allen.cells.rds'))  # Load allen cells data
  sce <- readRDS(get_path(variables$stored_variables, 'sce.preprocess.rds'))  # Load preprocessed SCE object
  mergeThem()  # Merge and split the data
  user.selection <- '4'  # Move to the next step
}

# Step 4: Feature selection and GLM correction
if (user.selection == '4') {
  cat(red("\nPerforming feature selection, adding chemistry information, feature exclusion and GLM correction\n"))
  cat("\n")
  sce <- readRDS(get_path(variables$stored_variables, 'sce.post.merge.rds'))  # Load merged SCE data
  feature_selection_and_correction()  # Perform feature selection and GLM correction
  user.selection <- '5'  # Move to the next step
}

# Step 5: Run the LIGER pipeline
if (user.selection == '5') {
  cat(red("\nRunning the liger pipeline\n"))
  cat("\n")
  sce <- readRDS(get_path(variables$stored_variables, 'preliger.rds'))  # Load data before LIGER
  features <- readRDS(get_path(variables$stored_variables, 'filtered.features.rds'))  # Load filtered features
  perform_liger()  # Execute LIGER analysis
  user.selection <- '6'  # Move to the next step
}

# Step 6: Execute nearest neighbor calculations
if (user.selection == '6') {
  cat(red("\nExecuting the nearest neighbor code\n"))
  cat("\n")
  sce <- readRDS(get_path(variables$stored_variables, 'sce.post.liger.rds'))  # Load post-LIGER SCE object
  calculate_nearest_neighbors()  # Compute nearest neighbors
  user.selection <- '7'  # Move to the next step
}

# Step 7: Filter clusters and exclude distant neighbors
if (user.selection == '7') {
  cat(red("\nPerforming cluster exclusion and distant neighbor exclusion\n"))
  cat("\n")
  sce <- readRDS(get_path(variables$stored_variables, 'sce.post.nearest.neighbors.rds'))  # Load SCE after nearest neighbors
  filter_sce()  # Filter out unwanted clusters or distant neighbors
  user.selection <- '8'  # Move to the next step
}

# Step 8: Perform SCmap layer 2/3 subtype assignment
if (user.selection == '8') {
  cat(red("\nPerforming layer 2/3 subtype assignment with scMap\n"))
  cat("\n")
  sce <- readRDS(get_path(variables$stored_variables, 'sce.filtered.rds'))  # Load filtered SCE object
  perform_scmap()  # Perform SCmap layer assignment
  user.selection <- '9'  # Move to the next step
}

# Step 9: Generate all related figures
if (user.selection == '9') {
  cat(red("\nGenerating the figures\n"))
  cat("\n")
  sce <- readRDS(get_path(variables$stored_variables, 'sce.filtered.rds'))  # Load filtered SCE object
  allen <- readRDS(get_path(variables$stored_variables, 'allen.cells.rds'))  # Load Allen cells data
  tenx.l23 <- readRDS(get_path(variables$stored_variables, 'l23.10x.post.scmap.rds'))  # Load 10x data post-SCmap
  get_the_figs()  # Generate all figures
  user.selection <- '10'  # Finish the analysis
}

# Print session information to file for record-keeping
cat("\n")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# Inform the user that the analysis is complete
cat(red("\nThe analysis has finished"))
cat("\n")
