Packages <- c('grid','Matrix','BiocParallel', 'scales', 'ggpubr','dplyr', 'extrafont','SingleCellExperiment',
             'Seurat','SeuratWrappers','pdftools','tidyverse','gplots','gridExtra','ggplot2','png','liger','tidyr')

suppressPackageStartupMessages(invisible(lapply(Packages, library , character.only = TRUE)));

library('ggplot2','Seurat','gridExtra','grid')
source('get_path.r')
source('single_cell_variables.r')

variables <- get_variables()

#get the layer 2/3 object

tenx.l23 <- readRDS(get_path(variables$stored_variables,'l23.10x.post.scmap.rds'))
tenx.l23.filtered <- subset(tenx.l23, subset = scmap_cluster_index == 'unassigned', invert = TRUE)

#asign metadata to layer 2/3 object
sample_vector <- tenx.l23.filtered[[]]$facs_sample
pc_type <- c()
stim_type <- c()
for (item in sample_vector){
    if ((grepl('1',item)) == 1) {
    pc_type <- append(pc_type,'low')
    }
    else if (grepl('2',item) == 1) {
    pc_type <- append(pc_type,'int')
    }
    else if (grepl('3',item) == 1) {
    pc_type <- append(pc_type,'int')
    }
    else{
    pc_type <- append(pc_type,'high')
    }
    if (grepl('vis',item) == 1) {
    stim_type <- append(stim_type,'vis')
    }
    else if (grepl('run',item) == 1) {
    stim_type <- append(stim_type,'run')
    }
    else {
    stim_type <- append(stim_type,'mm')
    }
}
tenx.l23.filtered$pc_type <- pc_type
tenx.l23.filtered$stim_type <- stim_type

# Set up the data
tenx.l23.run <- subset(tenx.l23.filtered, subset = stim_type == 'run')

tenx.l23.mm <- subset(tenx.l23.filtered, subset = stim_type == 'mm')
tenx.l23.mm.high <- subset(tenx.l23.mm, subset = pc_type == 'high')
tenx.l23.high.mm.umap.df <- as.data.frame(tenx.l23.mm.high[["umap"]]@cell.embeddings)

tenx.l23.vis <- subset(tenx.l23.filtered, subset = stim_type == 'vis')
tenx.l23.vis.high <- subset(tenx.l23.vis, subset = pc_type == 'high')
tenx.l23.vis.high.umap.df <- as.data.frame(tenx.l23.vis.high[["umap"]]@cell.embeddings)


adamts2_expressing <- subset(x = tenx.l23.filtered, subset = Adamts2 > 0);
tenx.l23.adamts2.umap.df <- as.data.frame(adamts2_expressing[["umap"]]@cell.embeddings)

baz1a_expressing <- subset(x = tenx.l23.filtered, subset = Baz1a > 0);
tenx.l23.baz1a.umap.df <- as.data.frame(baz1a_expressing[["umap"]]@cell.embeddings)

#Plotting function
get_l23_density_plot <- function(input_data,title,cVals,panel_label){
if (variables$FixAxes == TRUE){
    ggplot(input_data, aes(x=UMAP_1, y=UMAP_2)) +
    stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "Black",
                  bins = 10) +  # was 20
    scale_fill_distiller(palette = "Greys", direction = 1,limits = cVals,breaks = cVals, labels = cVals) +
    theme_classic() +
    ggtitle(title) +
    theme(text=element_text(family="ArialMT", size=12)) + theme(plot.title = element_text(hjust = 0.5)) +
    labs(y= "UMAP 2", x = "UMAP 1") + theme(legend.position="bottom") +
    scale_y_continuous(limits=c(-3,4), breaks = c(-3,4)) +
    scale_x_continuous(limits=c(-10,4), breaks = c(-10,4)) +
    theme(legend.position = c(0.15, 0.2), legend.direction = 'horizontal') +
    labs(fill = "Density") + labs(tag = panel_label) +
    theme(axis.ticks.length=unit(-0.25, "cm"))
}
else{
    ggplot(input_data, aes(x=UMAP_1, y=UMAP_2)) +
    stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "Black",
                  bins = 10) +  # was 20
    scale_fill_distiller(palette = "Greys", direction = 1,limits = cVals,breaks = cVals, labels = cVals) +
    theme_classic() +
    ggtitle(title) +
    theme(text=element_text(family="ArialMT", size=12)) + theme(plot.title = element_text(hjust = 0.5)) +
    labs(y= "UMAP 2", x = "UMAP 1") + theme(legend.position="right") +
    labs(fill = "Density") + labs(tag = panel_label) +
    theme(axis.ticks.length=unit(-0.25, "cm"))
    
}
        
}

#Get the photoconversion umaps
umap_pc_l23<- function(sce.input,title,panel_label){
Idents(sce.input) <- 'pc_type'
pc_colors <- c('low' = 'darkgreen','int' = 'gray', 'high' = 'red')
if (variables$FixAxes == TRUE){
DimPlot(sce.input, cols = pc_colors) + scale_y_continuous(limits=c(-3,4), breaks = c(-3,4)) +
    scale_x_continuous(limits=c(-10,4), breaks = c(-10,4)) + theme(legend.position = c(0.1, 0.9)) +
    theme(text=element_text(family="ArialMT", size=12)) + labs(y= "UMAP 2", x = "UMAP 1") +
    theme(axis.ticks.length=unit(-0.25, "cm")) + theme(axis.text = element_text(size = 12)) + ggtitle(title) +
     labs(tag = panel_label) 
}
else{
DimPlot(sce.input, cols = pc_colors) +
    theme(text=element_text(family="ArialMT", size=12)) + labs(y= "UMAP 2", x = "UMAP 1") +
    theme(axis.ticks.length=unit(-0.25, "cm")) + theme(axis.text = element_text(size = 12)) + ggtitle(title) +
     labs(tag = panel_label) 
}
}

#Get the plots
p1 <- umap_pc_l23(tenx.l23.mm,'Mismatch layer 2/3','A')
p2 <- umap_pc_l23(tenx.l23.vis,'Gratings layer 2/3','B')
p3 <- get_l23_density_plot(tenx.l23.high.mm.umap.df,'mismatch high',c(0, 0.1),'C')
p4 <- get_l23_density_plot(tenx.l23.vis.high.umap.df,'visual high',c(0, 0.1),'D')
p5 <- get_l23_density_plot(tenx.l23.adamts2.umap.df,'Adamts2 expressing',c(0, 0.1),'E')
p6 <- get_l23_density_plot(tenx.l23.baz1a.umap.df,'Baz1a expressing',c(0, 0.1),'F')

#Arrange and save them

panels = list(p1, p2, p3, p4, p5, p6)
             
lay <- rbind(c(1,1,1,1,2,2,2,2),
             c(1,1,1,1,2,2,2,2),
             c(1,1,1,1,2,2,2,2),
             c(3,3,3,3,4,4,4,4),
             c(3,3,3,3,4,4,4,4),
             c(3,3,3,3,4,4,4,4),
             c(5,5,5,5,6,6,6,6),
             c(5,5,5,5,6,6,6,6),
             c(5,5,5,5,6,6,6,6))

fig_path <- get_path(variables$figures,'fig_s2_c_through_h.pdf')
pdf(fig_path,width = 10, height = 15)
grid.arrange(grobs = panels, layout_matrix = lay)
dev.off()

grid.arrange(grobs = panels, layout_matrix = lay)