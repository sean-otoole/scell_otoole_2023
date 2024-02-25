get_the_figs <- function(){

allen.l23 <- subset(x = allen, idents = c('L2/3 IT'))

# color scheme formatting for umap plots
Idents(sce) <- 'categories'

sce.filtered.exclude.l23 <- subset(sce, subset = categories == 'L2/3 IT', invert = TRUE)
sce.filtered.exclude.l23$umap_labels <- sce.filtered.exclude.l23$categories
tenx.l23$umap_labels <- tenx.l23$scmap_cluster_index
sce.l23.sub.labels <- merge(tenx.l23, y = sce.filtered.exclude.l23, merge.dr = 'umap')

#setting umap colors

umap_colors <- c('L2/3 IT VISp Adamts2' = '#F24607', 'L2/3 IT VISp Agmat' = '#49D907',
                 'L2/3 IT VISp Rrad' = '#0597F2','L4' = 'grey','Lamp5' = 'grey','Macrophage' = 'grey',
                'Oligo' = 'grey', 'Pvalb' = 'grey','Sst' = 'grey','Vip' = 'grey',
                'VLMC' = 'grey')

#plotting individual panels
# global umap
Idents(sce.l23.sub.labels) <- 'umap_labels'
panel.a <- DimPlot(sce.l23.sub.labels, pt.size = 0.5, label = TRUE, cols = umap_colors) + theme(text=element_text(family="ArialMT"))+
xlab('Umap 1') + ylab('Umap 2') + NoLegend()

# adamts2/agmat/baz1a expression umaps
marker_positive <- WhichCells(object = tenx.l23, expression = Adamts2 > 0)
if (variables$FixAxes == TRUE){
panel.b <- DimPlot(tenx.l23, label=F, cells.highlight= marker_positive, cols.highlight = c("red"), cols= "grey")  +
scale_y_continuous(breaks = c(-2.5,0,3.5),limits = c(-2.5,3.5)) + 
scale_x_continuous(breaks = c(-6,0,2.0),limits = c(-6,2.0)) +
theme(text=element_text(family="ArialMT")) + NoLegend()
}
else{
panel.b <- DimPlot(tenx.l23, label=F, cells.highlight= marker_positive, cols.highlight = c("red"), cols= "grey")  +
theme(text=element_text(family="ArialMT")) + NoLegend()
}

marker_positive <- WhichCells(object = tenx.l23, expression = Agmat > 0)
    
if (variables$FixAxes == TRUE){
panel.c <- DimPlot(tenx.l23, label=F, cells.highlight= marker_positive, cols.highlight = c("green"), cols= "grey") +
scale_y_continuous(breaks = c(-2.5,0,3.5),limits = c(-2.5,3.5)) + 
scale_x_continuous(breaks = c(-6,0,2.0),limits = c(-6,2.0)) +
theme(text=element_text(family="ArialMT")) + NoLegend()
}
else{
panel.c <- DimPlot(tenx.l23, label=F, cells.highlight= marker_positive, cols.highlight = c("green"), cols= "grey") +
theme(text=element_text(family="ArialMT")) + NoLegend()
}
    
marker_positive <- WhichCells(object = tenx.l23, expression = Rrad > 0)
    
if (variables$FixAxes == TRUE){
panel.d <- DimPlot(tenx.l23, label=F, cells.highlight= marker_positive, cols.highlight = c("blue"), cols= "grey")  +
scale_y_continuous(breaks = c(-2.5,0,3.5),limits = c(-2.5,3.5)) + 
scale_x_continuous(breaks = c(-6,0,2.0),limits = c(-6,2.0)) +
theme(text=element_text(family="ArialMT")) + NoLegend()

}
else{
panel.d <- DimPlot(tenx.l23, label=F, cells.highlight= marker_positive, cols.highlight = c("blue"), cols= "grey")  +
theme(text=element_text(family="ArialMT")) + NoLegend()
}

###################### bootstrap code #########################

#simply gets the counts for the l23 cell types
format_from_seurat <- function(seurat_object,sample_column,cell_type_column){
    sample_names <- seurat_object[[sample_column]]
    cell_type_annotations <- seurat_object[[cell_type_column]]
    count_df <- cbind(sample_names,cell_type_annotations)
    not_Allen_index <- count_df[[sample_column]] != 'allen'
    count_df <- count_df[not_Allen_index,]
    count_df <- as.data.frame.matrix(table(count_df))
    count_df$l23_total <- rowSums(count_df[,c("L2/3 IT VISp Adamts2", "L2/3 IT VISp Agmat","L2/3 IT VISp Rrad")])
    if (variables$cosine_threshold == 0){
    count_df$total <- rowSums(count_df[1:(ncol(count_df))]) } else {
    count_df$total <- rowSums(count_df[1:(ncol(count_df)-1)]) #excludes the l23 total and counts unassigned
    }
    return(count_df)
}

#function to combine the middle conditions
combine_middle_conditions <- function(input_df){
input_df <- format_from_seurat(seurat_object,'facs_sample',current_mapping)
mm2a <- input_df['mm2a',]+input_df['mm3a',]
rownames(mm2a) <- 'mm2a'
mm2b <- input_df['mm2b',]+input_df['mm3b',]
rownames(mm2b) <- 'mm2b'
run2a<- input_df['run2a',]+input_df['run3a',]
rownames(run2a) <- 'run2a'
run2b <- input_df['run2b',]+input_df['run3b',]
rownames(run2b) <- 'run2b'
vis2a <- input_df['vis2a',]+input_df['vis3a',]
rownames(vis2a) <- 'vis2a'
vis2b <- input_df['vis2b',]+input_df['vis3b',]
rownames(vis2b) <- 'vis2b'
removal_list <- c('mm2a','mm3a','mm2b','mm3b','run2a','run3a','run2b','run3b','vis2a','vis3a',
                 'vis2b','vis3b')
input_df <- input_df[!(row.names(input_df) %in% removal_list),]
old_names <- rownames(input_df)
new_names <- gsub('4','3',old_names)
rownames(input_df) <- new_names
input_df <- rbind(input_df,mm2a,mm2b,run2a,run2b,vis2a,vis2b)
input_df <- input_df[ order(row.names(input_df)), ]
}

#for the simulation
get_average_l23 <- function(seurat_input){
    count_df <- format_from_seurat(seurat_input,'Sample',current_mapping)
    return(mean(count_df[,'l23_total']))
}


adjust_frequency <- function(input_df){
    input_df$adjustment_factor <- max(input_df$total)/input_df$total
    adjustment_factor <- input_df$adjustment_factor
    input_df <- input_df[,1:(ncol(input_df)-1)] * input_df$adjustment_factor
    input_df$adjustment_factor <- adjustment_factor
    return(input_df)
}

pool_after_adjustment <- function(input_df) {
    odd <- function(x) x%%2 != 0    #function for odd
    even <- function(x) x%%2 == 0   #function for even
    row_vector = c(1:nrow(input_df))    #generate a row vector
    a_set <- input_df[row_vector[odd(row_vector)],]    #take the odd values
    b_set <- input_df[row_vector[even(row_vector)],]   #take the even values
    added_set <- a_set[(1:nrow(a_set)-1),] + b_set    #df addition (with one row missing)
    missing_row <- a_set[nrow(a_set),] * 2   # compensate for only one vis sample set
    recombined_set <- rbind(added_set,missing_row)   #add the missing row
    old_names <- rownames(recombined_set)    #get the rownames
    new_names = substr(old_names,1,nchar(old_names)-1)   #ditch the last character
    rownames(recombined_set) <- new_names   #reassign
    return(recombined_set)
}


current_mapping <- 'scmap_cluster_index'
seurat_object <- tenx.l23
input_df <- format_from_seurat(seurat_object,'facs_sample',current_mapping)
input_df <- combine_middle_conditions(input_df)
input_df <- adjust_frequency(input_df)
input_df <- pool_after_adjustment(input_df)

# monte distribution functions

# esimating the fractions of l23 neurons from total sample
estimate_fractions <- function(input_df,condition){
    subset <- input_df[grep(condition, rownames(input_df)), ]
    Adam_sum <- sum(subset[,'L2/3 IT VISp Adamts2'])
    Agmat_sum <- sum(subset[,'L2/3 IT VISp Agmat'])
    Rrad_sum <- sum(subset[,'L2/3 IT VISp Rrad'])
    net_sum <- sum(Adam_sum, Agmat_sum, Rrad_sum)
    Adam_fraction <- Adam_sum/net_sum
    Agmat_fraction <- Agmat_sum/net_sum
    Rrad_fraction <- Rrad_sum/net_sum
    fractions <- c(Adam_fraction, Agmat_fraction, Rrad_fraction)
    names <- c('Adamts2','Agmat','Rrad')
    estimated_fractions <- rbind(names, fractions)
    return(estimated_fractions)
}

#generates names to choose from
generate_names_vector <- function(estimated_fractions,sample_size){
    adam_sample <- floor(sample_size * as.numeric(estimated_fractions[2,1]))
    agmat_sample <- floor(sample_size * as.numeric(estimated_fractions[2,2]))
    rrad_sample <- floor(sample_size * as.numeric(estimated_fractions[2,3]))
    adam_list <- rep('adam', times = adam_sample)
    agmat_list <- rep('agmat', times = agmat_sample)
    rrad_list <- rep('rrad', times = rrad_sample)
    full_list <- c(adam_list, agmat_list,rrad_list)
    return(full_list)
}

# gets the distribution of values from randomized data
generate_distrbution_dataframe <- function(names_vector,iteration,seurat_input){
    simulated_data_frame <- data.frame(adam=rep(0, iteration), agmat=rep(0, iteration), rrad=rep(0, iteration))
    simulated_sample_size <- get_average_l23(seurat_input)
    for (i in 1:iteration){
        set.seed(i)
        simulated_fractions <- proportions(table(sample(names_vector)[1:simulated_sample_size]))
        simulated_data_frame[i,] <- c(simulated_fractions['adam'],simulated_fractions['agmat'],simulated_fractions['rrad'])
    }
    return(simulated_data_frame)
}

# gets the fraction values
get_pc_fractions <- function(input_df,condition){
    subset_df <- input_df[grep(condition, rownames(input_df)), ]
    subset_df <- subset_df[,1:3]
    fraction_df <- subset_df/rowSums(subset_df)
    return(fraction_df)
}

# get the z values based on the simulations

get_z_df <- function(condition){

    z_df <- data.frame(Date=as.Date(character()),
                     File=character(), 
                     User=character(), 
                     stringsAsFactors=FALSE) 
    
    iterations <- variables$BootStrapIters
    total_fractions <-  estimate_fractions(input_df, condition)
    names_vector <- generate_names_vector(total_fractions,iterations)
    distribution_df <- generate_distrbution_dataframe(names_vector,iterations,seurat_object) #iter val is here
    condition_fractions <- get_pc_fractions(input_df, condition)
    cluster_list <- c('adam','agmat','rrad')
    fractions <- get_pc_fractions(input_df, condition)
    i = 0
    for (column in 1:(dim(fractions)[2])){
        i = i + 1
        current_cluster_fractions <- fractions[,column,drop = FALSE]
        current_z_df <- apply(current_cluster_fractions, 2, function(x) ((x - mean(distribution_df[,column]))/(sd(distribution_df[,column]))))
            if (i < 2){
                z_df <- current_z_df
            } else{
                z_df <- cbind(z_df, current_z_df)
            }
        
    }
    names(z_df)[1] <- paste('adam')
    names(z_df)[2] <- paste('agmat')
    names(z_df)[3] <- paste('rrad')
    return(z_df)
}

mm_z_df <- get_z_df('mm')
run_z_df <- get_z_df('run')
vis_z_df <- get_z_df('vis')

rrad_df <- cbind(mm_z_df[,3,drop = FALSE],run_z_df[,3,drop = FALSE],vis_z_df[,3,drop = FALSE]) 
colnames(rrad_df) <- c("Mismatch","Running","Visual")
agmat_df <- cbind(mm_z_df[,2,drop = FALSE],run_z_df[,2,drop = FALSE],vis_z_df[,2,drop = FALSE]) 
colnames(agmat_df) <- c("Mismatch","Running","Visual")
adam_df <- cbind(mm_z_df[,1,drop = FALSE],run_z_df[,1,drop = FALSE],vis_z_df[,1,drop = FALSE]) 
colnames(adam_df) <- c("Mismatch","Running","Visual")

# get the example distribution

example_dist <- function(seurat_object,iterations){
    input_df <- format_from_seurat(seurat_object,'Sample',current_mapping)   # originally second argument was 'Sample' which is just the paths
    input_df <- combine_middle_conditions(input_df)
    input_df <- adjust_frequency(input_df)
    input_df <- pool_after_adjustment(input_df)
    condition <- 'mm'
    i = 0
    total_fractions <-  estimate_fractions(input_df, condition)
    names_vector <- generate_names_vector(total_fractions,iterations)
    distribution_df <- generate_distrbution_dataframe(names_vector,iterations,seurat_object) #iter val is here
    adam <- (distribution_df$adam-mean(distribution_df$adam))/sd(distribution_df$adam)
    z_scores_df <- as.data.frame(adam)
    condition_fractions <- get_pc_fractions(input_df, condition)
    cluster <- 'adam'
    fractions <- get_pc_fractions(input_df, condition)
    current_column_number <- grep(pattern=cluster, x=colnames(fractions), ignore.case = TRUE)
    current_column <- fractions[,current_column_number]
    current_column <- (current_column-mean(distribution_df$adam))/sd(distribution_df$adam)
    example_dist_fig <- ggplot(z_scores_df, aes_string(x='adam')) + geom_density(color="black", fill="gray") +
    geom_vline(xintercept=current_column[1], col = "gray") +
    geom_vline(xintercept=current_column[2], col = "pink") +
    geom_vline(xintercept=current_column[3], col = "red") +
    geom_vline(xintercept=2.58, col = "gray", size = 1, linetype = 'dotted') +   #pval of 0.01
    geom_vline(xintercept=-2.58, col = "gray", size = 1,linetype = 'dotted') +
    ggtitle(paste(condition,cluster, sep=" ", collapse=NULL)) +
    theme_classic() +
    labs(y = "Density", x = "Fraction of L2/3 Neurons") +
    theme(axis.text.y=element_blank()) +
    theme(text=element_text(family="ArialMT")) +
    ggtitle("") + xlim(-10, 10) +
    theme(text=element_text(size=30))
    return(example_dist_fig)
}

mm_dist <- example_dist(tenx.l23, variables$BootStrapIters)

# get the monte dotplots
shiny_monte_plot <- function(values){
pc = rep(c('PC1','PC2','PC3'),3)
clusters = c(rep('Mismatch',3),rep('Running',3),rep("Visual",3))
test_df <- data.frame(z_value = values, photoconversion = pc, clusters = clusters)
monte_plot <- ggplot(test_df, aes(x = clusters, y=z_value,color = photoconversion)) +  geom_point(size = 5) +
ylab('Z score') + theme(axis.title.y = element_text(size=30)) +
theme(axis.title.x = element_text(size=30)) +
scale_color_manual(breaks = c("PC1", "PC2", "PC3"),values=c("gray", "pink","red")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
geom_hline(yintercept=0, col = "black", size = 1) + ylim(-10, 10) +
geom_hline(yintercept=2.58, col = "gray", size = 1, linetype = 'dotted') +   #pval of 0.01
geom_hline(yintercept=-2.58, col = "gray", size = 1,linetype = 'dotted') +
theme(axis.text.x = element_text(angle = 70, vjust = 0.97, hjust=1.0)) +
theme(axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm")) +
theme(axis.title.x=element_blank())  +
theme(legend.position="right") +
theme(axis.line = element_line(colour = 'black', size = 1)) +
theme(axis.ticks = element_line(colour = "black", size = 1)) +
theme(legend.position = "none")
return(monte_plot)
}

mismatch_monte <- shiny_monte_plot(c((mm_z_df[,1]),(mm_z_df[,2]),(mm_z_df[,3])))
run_monte <- shiny_monte_plot(c((run_z_df[,1]),(run_z_df[,2]),(run_z_df[,3]))) 
vis_monte <- shiny_monte_plot(c((vis_z_df[,1]),(vis_z_df[,2]),(vis_z_df[,3])))

adam_monte <- shiny_monte_plot(c((adam_df[,1]),(adam_df[,2]),(adam_df[,3]))) + ggtitle('Adamts2') +
theme(text=element_text(size=30)) + theme(plot.title = element_text(size=30)) +
theme(plot.title = element_text(hjust = 0.5)) +
theme(text=element_text(family="ArialMT"))

agmat_monte <- shiny_monte_plot(c((agmat_df[,1]),(agmat_df[,2]),(agmat_df[,3])))+ ggtitle('Agmat') +
theme(text=element_text(size=30)) + theme(plot.title = element_text(size=30)) +
theme(plot.title = element_text(hjust = 0.5)) +
theme(text=element_text(family="ArialMT"))

rrad_monte <- shiny_monte_plot(c((rrad_df[,1]),(rrad_df[,2]),(rrad_df[,3])))+ ggtitle('Rrad') +
theme(text=element_text(size=30)) + theme(plot.title = element_text(size=30)) +
theme(plot.title = element_text(hjust = 0.5)) +
theme(text=element_text(family="ArialMT"))


panel.e <- mm_dist
panel.f <- adam_monte
panel.g <- agmat_monte
panel.h <- rrad_monte

###################### figure 2 and 4 assembly code #####################################

panels = list(panel.a,panel.b,panel.c,panel.d,panel.e,panel.f,panel.g,panel.h)
             
lay <- rbind(c(1,1,1,1,1,1,1,1,2,2,2,2),
             c(1,1,1,1,1,1,1,1,2,2,2,2),
	     c(1,1,1,1,1,1,1,1,3,3,3,3),
	     c(1,1,1,1,1,1,1,1,3,3,3,3),
	     c(1,1,1,1,1,1,1,1,4,4,4,4),
	     c(1,1,1,1,1,1,1,1,4,4,4,4),
	     c(5,5,5,6,6,6,7,7,7,8,8,8),
	     c(5,5,5,6,6,6,7,7,7,8,8,8),
             c(5,5,5,6,6,6,7,7,7,8,8,8))

fig3_path <- get_path(variables$figures,"figure2_and_4.pdf")
pdf(fig3_path, width = 20, height = 20)
grid.arrange(grobs = panels, layout_matrix = lay)
dev.off()

################## figure 3 A and B and supplements code ##############################

Idents(tenx.l23) <- 'scmap_cluster_index'

if (nrow(unique(tenx.l23[[]]$scmap_cluster_index)) == 3){
	cat(red("\nNo unassiagned l23 types\n")) } else {
	tenx.l23 <- subset(x = tenx.l23, idents = c('unassigned'), invert = TRUE)
}

# table function for barplot panels

format_for_stacked_bar <- function(sce_object,sample_identities,cell_type_identities){
samples <- sce_object[[]] %>% pull(!!sym(sample_identities))
categories <- sce_object[[]] %>% pull(!!sym(cell_type_identities))
df <- data.frame(samples, categories)
count.df <- count(df, samples, categories)
tmp <- complete(count.df,samples, categories) # fill in the missing values
tmp <-as.data.frame(tmp)  # ditch the tibble format
iteration_list <- 1:nrow(tmp)
tmp[is.na(tmp)] <- 0  # convert nas to zeroes
percentage.values <- lapply(iteration_list, function(x) ((tmp[x,'n'])/sum(tmp[(tmp$samples == tmp[x,'samples']),]$n))*100)
tmp$percentages <- unlist(percentage.values)
tmp$n <- NULL
sample_list <- tmp$samples
sample_list <- lapply(sample_list, function(x) substr(x,1,nchar(x)-1))
sample_list <- str_replace(sample_list,'vis3','vis2')
sample_list <- str_replace(sample_list,'vis4','vis3')
sample_list <- str_replace(sample_list,'run3','run2')
sample_list <- str_replace(sample_list,'run4','run3')
sample_list <- str_replace(sample_list,'mm3','mm2')
sample_list <- str_replace(sample_list,'mm4','mm3')       
tmp$samples <- sample_list
iteration_list <- 1:nrow(tmp)
combined_names <- lapply(iteration_list, function(x) paste(tmp$samples[x],tmp$categories[x],sep = "-"))
tmp$combined <- unlist(combined_names)
tmp <- tmp %>% group_by(combined) %>% summarise(percentage = mean(percentages))
tmp <- as.data.frame(tmp)
tmp <- tmp %>% separate(combined, c("samples", "categories"), sep = "-")
return(tmp)
}

# dot plot, panel a
Idents(object = sce) <- 'categories'
levels(sce) <- rev(c('Macrophage','Oligo','L4','L2/3 IT','Sst','Pvalb','Vip','Lamp5'))
markers.to.plot <- c('Cd68','Olig1','Sst','Pvalb','Vip','Rspo1','Rorb','Lamp5','Wfs1','Gad2','Gad1','Cux2','Slc17a7')
panel.a <- DotPlot(sce, features = rev(markers.to.plot), dot.scale = 10, scale.max = 50,
                  col.min = 0, col.max = 1) + RotatedAxis() + theme(text=element_text(family="ArialMT")) + theme(legend.position= "right")+
		  ggtitle(sprintf("tenx cell count:%d", nrow(sce[[]])))

# global types distribution, panel b
tmp <- format_for_stacked_bar(sce,'facs_sample','categories')
nb.cols <- length(unique(tmp$categories))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlGn"))(nb.cols)
panel.b <- ggplot(tmp, aes(x = samples, y = percentage, fill = categories)) + geom_bar(position = "fill",stat = "identity") +
scale_y_continuous(labels = scales::percent_format()) +
scale_fill_manual(values = mycolors) +
theme_classic() +
RotatedAxis()


#l23 subtype distribution, panel c
tmp <- format_for_stacked_bar(tenx.l23,'facs_sample','umap_labels')
nb.cols <- length(unique(tmp$categories))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlGn"))(nb.cols)
panel.c <- ggplot(tmp, aes(x = samples, y = percentage, fill = categories)) + geom_bar(position = "fill",stat = "identity") +
scale_y_continuous(labels = scales::percent_format()) +
scale_fill_manual(values = mycolors) +
theme_classic() +
RotatedAxis()

# l23 heatmap plot, panels D and E
heat_map_markers <- c('Gsta4','Timp2','Pea15a','Rfx3','Nov','Matn2','Cdh13','Adamts2',     
                      'Slc24a2','Scn1b','Rasgrp1','AI593442','Cdh12','Igfn1','AI593442','Agmat',
                      'Grasp','Nptx2','Egr4','Tnfaip6','Mas1','Ptgs2','Nefl','Rrad','Baz1a')
heat.min = -1
heat.max = 1
tenx.l23 <- ScaleData(object = tenx.l23, features = heat_map_markers,vars.to.regress = c("nCount_RNA","percent.mt","facs_sample","chemistry","nFeature_RNA"))
Idents(tenx.l23) <- 'scmap_cluster_index'
adam <- subset(x = tenx.l23, idents = c('L2/3 IT VISp Adamts2'))
agmat <- subset(x = tenx.l23, idents = c('L2/3 IT VISp Agmat'))
rrad <- subset(x = tenx.l23, idents = c('L2/3 IT VISp Rrad'))
scdm.adam <- GetAssayData(object = adam, slot = 'scale.data')
adam.names <- colnames(scdm.adam)
scdm.agmat <- GetAssayData(object = agmat, slot = 'scale.data')
agmat.names <- colnames(scdm.agmat)
scdm.rrad <- GetAssayData(object = rrad, slot = 'scale.data')
rrad.names <- colnames(scdm.rrad)
rolling_window <- 20
scdm.adam <- t(rollmean(t(scdm.adam),rolling_window, na.pad = TRUE, partial = TRUE, fill = NA))
scdm.agmat <- t(rollmean(t(scdm.agmat),rolling_window, na.pad = TRUE, partial = TRUE, fill = NA))
scdm.rrad <- t(rollmean(t(scdm.rrad),rolling_window, na.pad = TRUE, partial = TRUE, fill = NA))
scdm.all <- cbind(scdm.adam,scdm.agmat,scdm.rrad)
colnames(scdm.all) <- c(adam.names,agmat.names,rrad.names)
tenx.l23.temp <- tenx.l23
tenx.l23.temp <- SetAssayData(tenx.l23.temp, slot = 'scale.data', scdm.all)
panel.d <- DoHeatmap(tenx.l23.temp, features = heat_map_markers, label = FALSE,disp.min=heat.min, disp.max = heat.max,
         raster = FALSE, draw.lines = FALSE) +
scale_fill_gradientn(colours = c("blue", "white", "red"),na.value = "white",values = scales::rescale(c(-3.0, -1.0, 0, 1.0, 3.0)))+
ggtitle(sprintf("tenx cell count:%d", nrow(tenx.l23.temp[[]])))

allen.l23 <- ScaleData(object = allen.l23, features = heat_map_markers,vars.to.regress = c("nCount_RNA","percent.mt"))
Idents(object = allen.l23) <- 'allenCluster'
levels(allen.l23) <- c('L2/3 IT VISp Adamts2','L2/3 IT VISp Agmat','L2/3 IT VISp Rrad')
panel.e <- DoHeatmap(allen.l23, features = heat_map_markers, label = FALSE,disp.min=-2.5, disp.max = 2.5,
         raster = FALSE,draw.lines = FALSE) +
scale_fill_gradientn(colours = c("blue", "white", "red"),values = scales::rescale(c(-3.0, -1.0, 0, 1.0, 3.0)))+
ggtitle(sprintf("allen cell count:%d", nrow(allen.l23[[]])))

#figure generation code
panels = list(panel.a,panel.b,panel.c,panel.d,panel.e)             
lay <- rbind(c(1,1,1,1,1,1,2,2,2,3,3,3),
             c(1,1,1,1,1,1,2,2,2,3,3,3),
             c(1,1,1,1,1,1,2,2,2,3,3,3),
             c(1,1,1,1,1,1,2,2,2,3,3,3),
             c(4,4,4,4,4,4,5,5,5,5,5,5),
             c(4,4,4,4,4,4,5,5,5,5,5,5),
             c(4,4,4,4,4,4,5,5,5,5,5,5),
             c(4,4,4,4,4,4,5,5,5,5,5,5))

fig2_path <- get_path(variables$figures,"figure_3_AB_and_supplements.pdf")
pdf(fig2_path, width = 20, height = 20)
grid.arrange(grobs = panels, layout_matrix = lay)
dev.off()

################## summary stats code ##############################

sample_list <- c('vis1a','vis1b','vis2a','vis2b','vis3a','vis3b','vis4b',
                'mm1a','mm1b','mm2a','mm2b','mm3a','mm3b','mm4a','mm4b',
                'run1a','run1b','run2a','run2b','run3a','run3b','run4a','run4b')

base_directory <- '/home_fmi/01/otoosean/home_faimsrv01/Sequel/big_project/'
path_add_on <- '/outs/metrics_summary.csv'

build_summary_path <-function(condition){
    paste(base_directory,condition,path_add_on, sep = '')
}

sample_summary_paths <- lapply(sample_list,build_summary_path)

sampleSummary <- data.frame()

for (path in sample_summary_paths){
    current_summary <- read.csv(path)
    rownames(current_summary) <- qdapRegex::ex_between(path, base_directory, path_add_on)[[1]]
    sampleSummary <- rbind(sampleSummary,current_summary)
}

sampleSummary[] <- lapply(sampleSummary, function(x) as.numeric(gsub("[$,%]", "", x)))

sampleSummary$Chemistry <- c('v2','v3','v2','v3','v2','v3','v3',
                             'v3','v3','v3','v3','v3','v3','v3','v3',
                            'v2','v3','v2','v3','v2','v3','v2','v3')


x.max <- max(sampleSummary$Number.of.Reads)
y.max <- max(sampleSummary$Median.Genes.per.Cell)

features_vs_readCount <- ggscatter(sampleSummary, x = "Number.of.Reads", 
                                   y = "Median.Genes.per.Cell",
                                   color = "Chemistry", 
                                   palette = "jco", xlim = c(0,x.max),
                                   ylim = c(0,y.max), size = 3, xlab = 'Read Count', 
                                   ylab = 'Feature Count', show.legend.text = TRUE) +
                                   theme(axis.ticks.length=unit(-0.1, "cm"))

depth_vs_saturation <- ggscatter(sampleSummary, x = "Number.of.Reads", 
                                   y = "Sequencing.Saturation",
                                   color = "Chemistry", 
                                   ylim = c(0,100), xlim = c(0,max(sampleSummary$Number.of.Reads)),
                                   palette = "jco", size = 3, xlab = 'Read depth', 
                                   ylab = 'Sequencing saturation', show.legend.text = TRUE) +
                                   theme(axis.ticks.length=unit(-0.1, "cm"))

features_vs_umiCount <- ggscatter(sampleSummary, x = 'Median.Genes.per.Cell', 
                                   y = "Mean.Reads.per.Cell",
                                   color = "Chemistry", 
                                         palette = "jco",
                                   xlim = c(0,max(sampleSummary$Median.Genes.per.Cell)), 
                                   ylim = c(0,max(sampleSummary$Mean.Reads.per.Cell)),
                                   size = 3, xlab = 'median genes per cell', 
                                   ylab = 'Median UMI Count', show.legend.text = TRUE) +
                                   theme(axis.ticks.length=unit(-0.1, "cm"))


readCount_vs_genes.detected <- ggscatter(sampleSummary, x = "Number.of.Reads", 
                                         y = "Total.Genes.Detected",
                                         color = "Chemistry", 
                                         palette = "jco", xlim = c(0,max(sampleSummary$Number.of.Reads)),
                                         ylim = c(16000,max(sampleSummary$Total.Genes.Detected)),
                                         size = 3, xlab = 'Read depth', 
                                         ylab = 'Genes Detected', show.legend.text = TRUE) +
                                         theme(axis.ticks.length=unit(-0.1, "cm"))

                                   
features_vs_readCount <- arrangeGrob(features_vs_readCount, top = textGrob("A", x = unit(0, "npc")
         , y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18, fontfamily="ArialMT")))

depth_vs_saturation <- arrangeGrob(depth_vs_saturation, top = textGrob("B", x = unit(0, "npc")
         , y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18, fontfamily="ArialMT")))

features_vs_umiCount <- arrangeGrob(features_vs_umiCount, top = textGrob("D", x = unit(0, "npc")
         , y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18, fontfamily="ArialMT")))

readCount_vs_genes.detected <- arrangeGrob(readCount_vs_genes.detected, top = textGrob("C", x = unit(0, "npc")
         , y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18, fontfamily="ArialMT")))

g <- grid.arrange(features_vs_readCount, depth_vs_saturation,readCount_vs_genes.detected,features_vs_umiCount, nrow = 2, ncol=2)

scell_summary_fig_path <- get_path(variables$figures,"single_cell_summary_stats.pdf")
ggsave(scell_summary_fig_path , plot = g, device=cairo_pdf, width=20, height=20)

################# coverage plots  #############################################

mm1.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('mm1a','mm1b')) 
mm1.plt <- DimPlot(mm1.cov) + ggtitle('mm low')
mm2.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('mm2a','mm2b','mm3a','mm3b'))
mm2.plt <- DimPlot(mm2.cov) + ggtitle('mm low') + ggtitle('mm mid')
mm3.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('mm4a','mm4b'))
mm3.plt <- DimPlot(mm3.cov) + ggtitle('mm low') + ggtitle('mm high')

run1.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('run1a','run1b'))
run1.plt <- DimPlot(run1.cov) + ggtitle('run low')
run2.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('run2a','run2b','run3a','run3b'))
run2.plt <- DimPlot(run2.cov) + ggtitle('run mid')
run3.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('run4a','run4b'))
run3.plt <- DimPlot(run3.cov) + ggtitle('run high')

vis1.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('vis1a','vis1b'))
vis1.plt <- DimPlot(vis1.cov) + ggtitle('vis low')
vis2.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('vis2a','vis2b','vis3a','vis3b'))
vis2.plt <- DimPlot(vis2.cov) + ggtitle('vis mid')
vis3.cov <- subset(sce.l23.sub.labels, subset = facs_sample %in% c('vis4b'))
vis3.plt <- DimPlot(vis3.cov) + ggtitle('vis high')

panels = list(mm1.plt,mm2.plt,mm3.plt,run1.plt,run2.plt,run3.plt,vis1.plt,vis2.plt,vis3.plt)             
lay <- rbind(c(1,1,1,2,2,2,3,3,3),
             c(1,1,1,2,2,2,3,3,3),
             c(1,1,1,2,2,2,3,3,3),
             c(4,4,4,5,5,5,6,6,6),
             c(4,4,4,5,5,5,6,6,6),
             c(4,4,4,5,5,5,6,6,6),
             c(7,7,7,8,8,8,9,9,9),
             c(7,7,7,8,8,8,9,9,9),
             c(7,7,7,8,8,8,9,9,9))

coverage_plots_path <- get_path(variables$figures,"single_cell_coverage_plots.pdf")
pdf(coverage_plots_path , width = 20, height = 20)
grid.arrange(grobs = panels, layout_matrix = lay)
dev.off()

################# p_values and l23 sample sizes  ##############################

library("gridExtra")
monte_p_values <- get_path(variables$figures,"monte_p_values.pdf")
pdf(monte_p_values)
grid.table(rbind(signif(2*pnorm(-abs(mm_z_df)),2),
      signif(2*pnorm(-abs(run_z_df)),2),
      signif(2*pnorm(-abs(vis_z_df)),2)))
dev.off()

sample_list <- tenx.l23$facs_sample
cell_types <- tenx.l23$scmap_cluster_index
sample_list <- lapply(sample_list, function(x) substr(x,1,nchar(x)-1))
sample_list <- str_replace(sample_list,'vis3','vis2')
sample_list <- str_replace(sample_list,'vis4','vis3')
sample_list <- str_replace(sample_list,'run3','run2')
sample_list <- str_replace(sample_list,'run4','run3')
sample_list <- str_replace(sample_list,'mm3','mm2')
sample_list <- str_replace(sample_list,'mm4','mm3') 
sample.celltypes.df <- data.frame(sample_list,cell_types)
sample.celltypes.df <- as.data.frame.matrix(table(sample.celltypes.df))
sample.celltypes.df[nrow(sample.celltypes.df)+1,] <- colSums(sample.celltypes.df[,])
row.names(sample.celltypes.df)[10] <- "sum" 
l23_samples_sizes_pdf <- get_path(variables$figures,"sample_sizes_layer_23.pdf")
pdf(l23_samples_sizes_pdf)
grid.table(sample.celltypes.df)
dev.off()

# report the l23 cell fractions
fractions_df <- format_for_stacked_bar(tenx.l23,'facs_sample','umap_labels')
colnames(fractions_df) <- c('sample','cell type','fraction')
fractions_df$fraction <- fractions_df$fraction/100

fractions_df_path <- get_path(variables$figures,"fractions_df.rds")
saveRDS(fractions_df, fractions_df_path)
fractions_df[,3] <-round(fractions_df[,3],4)

l23_fractions_pdf <- get_path(variables$figures,"fractions_layer_23.pdf")
pdf(l23_fractions_pdf, height = 10) 
grid.table(fractions_df)
dev.off()

################# variables used  #############################################

variables_summary_path = get_path(variables$figures,'variables_used.pdf')
filtered_features_path = get_path(variables$stored_variables,'filtered.features.rds')

pdf(variables_summary_path)
plot.new()
text(x=.1, y=.9, adj = 0,sprintf("mito_max:%f", variables$mito_max))
text(x=.1, y=.85, adj = 0,sprintf("mito_min:%f", variables$mito_min))
text(x=.1, y=.8, adj = 0,sprintf("FDR_thresh:%f", variables$FDR_thresh))
text(x=.1, y=.75, adj = 0,sprintf("cosine_threshold:%f", variables$cosine_threshold))
text(x=.1, y=.7, adj = 0,sprintf("initial_features:%f", variables$initial_features))
text(x=.1, y=.65, adj = 0,sprintf("k_value:%f", variables$k_value))
text(x=.1, y=.6, adj = 0,sprintf("final_features:%f", length(readRDS(filtered_features_path ))))
text(x=.1, y=.55, adj = 0,sprintf("exclude_features:%s", variables$exclude_features))
text(x=.1, y=.5, adj = 0,sprintf("lambda_value:%f", variables$lambda_value))
text(x=.1, y=.45, adj = 0,sprintf("nearest_neighbors_cutoff:%f", variables$nearest_neighbors_cutoff))
text(x=.1, y=.4, adj = 0,sprintf("removed_dims:%f", variables$removed_dims))
text(x=.1, y=.35, adj = 0,sprintf("neighbor_values:%f", variables$neighbor_values))
dev.off()
################################################################################

### generate the remaining figures from additional scripts

source('fig_s2b.r')
cat(red("\nFinished generating figure S2B\n"))
cat("\n")

source('fig_s3_b_through_e.r')
cat(red("\nFinished generating figure S3B through E\n"))
cat("\n")

source('fig_s2_c_through_h.r')
cat(red("\nFinished generating figure S2C through H\n"))
cat("\n")

source('fig_s4L.r')
cat(red("\nFinished generating figure S4L\n"))
cat("\n")

source('fig_s5.r')
cat(red("\nFinished generating figure S5L\n"))
cat("\n")

source('fig_3c.r')
cat(red("\nFinished generating figure 3CL\n"))
cat("\n")

}