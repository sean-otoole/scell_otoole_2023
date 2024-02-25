library('ggplot2')
library('Seurat')
library('dplyr')
library('extrafont')

source('get_path.r')
source('single_cell_variables.r')

variables <- get_variables()

allen.cells <- readRDS('stored_variables/allen.cells.rds')

#get allen plot
allen.l23 <- subset(x = allen.cells, idents = c('L2/3 IT'))
Idents(allen.l23) <- 'allenCluster'
fig_path <- get_path(variables$figures,'fig_s2b.pdf')

if (variables$FixAxes == TRUE){
allen_plot <- DimPlot(allen.l23) + 
scale_y_continuous(breaks = c(-3,0,3),limits = c(-3,3)) + 
scale_x_continuous(breaks = c(9,15),limits = c(9,15)) +
theme(text=element_text(family="ArialMT")) + NoLegend()
} else {
allen_plot <- DimPlot(allen.l23, label = TRUE) +
theme(text=element_text(family="ArialMT")) + NoLegend()
}

ggsave(fig_path,allen_plot)