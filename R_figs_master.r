suppressPackageStartupMessages((library('crayon')))

cat(red("\nLoading packages this will take a few seconds\n"))

Packages <- c("DropletUtils", "scater", "TxDb.Mmusculus.UCSC.mm10.ensGene",'lsa','zoo','grid',
             'AnnotationHub','scRNAseq','scran','Matrix','BiocParallel', 'scales', 'ggpubr',
             'BiocSingular','dplyr','edgeR','pheatmap', 'extrafont','SingleCellExperiment',
             'Seurat','stringr','SeuratWrappers','pdftools','tidyverse','gplots','gridExtra',
	     'dendextend','RColorBrewer','ggplot2','png','liger','rliger','scmap','tidyr','stargazer','qdapRegex')

suppressPackageStartupMessages(invisible(lapply(Packages, library , character.only = TRUE)));

cat(red("\nFinished with package loading\n"))

source('get_path.r')
source('single_cell_variables.r')
source('generate_allen_data.r')
source('preprocess.r')
source('merge_and_split.r')
source('feature_selection_and_correction.r')
source('perform_liger.r')
source('calculate_nearest_neighbors.r')
source('filter_sce.r')
source('l23_scmap_assignment.r')
source('get_the_figs.r')

variables <- get_variables()

#sample_list <- c('vis3a') # just for testing

sample_list <- c('vis1a','vis1b','vis2a','vis2b','vis3a','vis3b','vis4b',
                'mm1a','mm1b','mm2a','mm2b','mm3a','mm3b','mm4a','mm4b',
                'run1a','run1b','run2a','run2b','run3a','run3b','run4a','run4b')

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

cat("\n")
cat("\n")

user.selection <- readLines(con="stdin", 1)

#cat(fil, "\n")

if (user.selection == '1') {
cat(red("\nGenerating the allen data"))
cat("\n")
generate_allen_data()
user.selection <- '2'
}

if (user.selection == '2') {
cat(red("\nPre-processing the data"))
cat("\n")
preprocess_pipeline()
user.selection <- '3'
} 

if (user.selection == '3') {
cat(red("\nPerforming merge and split\n"))
cat("\n")
allen.cells <- readRDS(get_path(variables$stored_variables,'allen.cells.rds'))
sce <- readRDS(get_path(variables$stored_variables,'sce.preprocess.rds'))
mergeThem()
user.selection <- '4'
}

if (user.selection == '4') {
cat(red("\nPerforming feature selection, adding chemistry information,feature exclusion and GLM correction\n"))
cat("\n")
sce <- readRDS(get_path(variables$stored_variables,'sce.post.merge.rds'))
feature_selection_and_correction()
user.selection <- '5'
}

if (user.selection == '5') {
cat(red("\nRunning the liger pipeline\n"))
cat("\n")
sce <- readRDS(get_path(variables$stored_variables,'preliger.rds'))
features <- readRDS(get_path(variables$stored_variables,'filtered.features.rds'))
perform_liger()
user.selection <- '6'
}

if (user.selection == '6') {
cat(red("\nExecuting the nearest neighbor code\n"))
cat("\n")
sce <- readRDS(get_path(variables$stored_variables,'sce.post.liger.rds'))
calculate_nearest_neighbors()
user.selection <- '7'
}

if (user.selection == '7') {
cat(red("\nPerforming cluster exclusion and distant neighbor exclusion\n"))
cat("\n")
sce <- readRDS(get_path(variables$stored_variables,'sce.post.nearest.neighbors.rds'))
filter_sce()
user.selection <- '8'
}

if (user.selection == '8') {
cat(red("\nPerforming layer 2/3 subtype assignment with scMap\n"))
cat("\n")
sce <- readRDS(get_path(variables$stored_variables,'sce.filtered.rds'))
perform_scmap()
user.selection <- '9'
}

if (user.selection == '9') {
cat(red("\nGenerating the figures\n"))
cat("\n")
sce <- readRDS(get_path(variables$stored_variables,'sce.filtered.rds'))
allen <- readRDS(get_path(variables$stored_variables,'allen.cells.rds'))
tenx.l23 <- readRDS(get_path(variables$stored_variables,'l23.10x.post.scmap.rds'))
get_the_figs()
user.selection <- '10'
}

cat("\n")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat(red("\nThe analysis has finished"))

cat("\n")

