Packages <- c("DropletUtils", "scater", "TxDb.Mmusculus.UCSC.mm10.ensGene",'lsa','zoo','grid',
             'AnnotationHub','scRNAseq','scran','Matrix','BiocParallel', 'scales', 'ggpubr',
             'BiocSingular','dplyr','edgeR','pheatmap', 'extrafont','SingleCellExperiment',
             'Seurat','stringr','SeuratWrappers','pdftools','tidyverse','gplots','gridExtra',
	     'dendextend','RColorBrewer','ggplot2','png','liger','rliger','scmap','tidyr','stargazer')

suppressPackageStartupMessages(invisible(lapply(Packages, library , character.only = TRUE)));
source('get_path.r')
source('single_cell_variables.r')

variables <- get_variables()

sce <- readRDS(get_path(variables$stored_variables,'sce.filtered.rds'))

# kick out everyone except major interneuon groups
exclusion_list <- c('Astro','Meis2','Peri','Serpinf1','SMC','Sncg','VLMC','L2/3 IT',
                   'L4','Macrophage','Oligo')
sce <- subset(sce, subset = categories %in% exclusion_list, invert = TRUE)


###################### monte code #########################

#simply gets the counts for the cell groups
format_from_seurat <- function(seurat_object,sample_column,cell_type_column){
    sample_names <- seurat_object[[sample_column]]
    cell_type_annotations <- seurat_object[[cell_type_column]]
    count_df <- cbind(sample_names,cell_type_annotations)
    not_Allen_index <- count_df[[sample_column]] != 'allen'
    count_df <- count_df[not_Allen_index,]
    count_df <- as.data.frame.matrix(table(count_df))
    count_df$group_total <- rowSums(count_df[,c("Lamp5","Pvalb","Sst","Vip")])
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
get_average_group <- function(seurat_input){
    count_df <- format_from_seurat(seurat_input,'Sample',current_mapping)
    return(mean(count_df[,'group_total']))
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


current_mapping <- 'categories'
seurat_object <- sce
input_df <- format_from_seurat(seurat_object,'facs_sample',current_mapping)
input_df <- combine_middle_conditions(input_df)
input_df <- adjust_frequency(input_df)
input_df <- pool_after_adjustment(input_df)

# bootstrap distribution functions

# esimating the fractions of l23 neurons from total sample
estimate_fractions <- function(input_df,condition){
    subset <- input_df[grep(condition, rownames(input_df)), ]
    lamp5_sum <- sum(subset[,'Lamp5'])
    pvalb_sum <- sum(subset[,'Pvalb'])
    sst_sum <- sum(subset[,'Sst'])
    vip_sum <- sum(subset[,'Vip'])
    net_sum <- sum(lamp5_sum,pvalb_sum,sst_sum,vip_sum)
    lamp5_fraction <- lamp5_sum/net_sum
    pvalb_fraction <- pvalb_sum/net_sum
    sst_fraction <- sst_sum/net_sum
    vip_fraction <- vip_sum/net_sum
    fractions <- c(lamp5_fraction,pvalb_fraction,sst_fraction,vip_fraction)
    names <- c('lamp5','pvalb','sst','vip')
    estimated_fractions <- rbind(names, fractions)
    return(estimated_fractions)
}

#generates names to choose from
generate_names_vector <- function(estimated_fractions,sample_size){
    lamp5_sample <- floor(sample_size * as.numeric(estimated_fractions[2,1]))
    pvalb_sample <- floor(sample_size * as.numeric(estimated_fractions[2,2]))
    sst_sample <- floor(sample_size * as.numeric(estimated_fractions[2,3]))
    vip_sample <- floor(sample_size * as.numeric(estimated_fractions[2,4]))
    lamp5_list <- rep('lamp5', times = lamp5_sample)
    pvalb_list <- rep('pvalb', times = pvalb_sample)
    sst_list <- rep('sst', times = sst_sample)
    vip_list <- rep('vip', times = vip_sample)
    full_list <- c(lamp5_list,pvalb_list,sst_list,vip_list)
    return(full_list)
}

# gets the distribution of values from randomized data
generate_distrbution_dataframe <- function(names_vector,iteration,seurat_input){
    simulated_data_frame <- data.frame( lamp5=rep(0, iteration),pvalb=rep(0, iteration),
                                       sst=rep(0, iteration),vip=rep(0, iteration))
    simulated_sample_size <- get_average_group(seurat_input)
    for (i in 1:iteration){
        set.seed(i)
        simulated_fractions <- proportions(table(sample(names_vector)[1:simulated_sample_size]))
        simulated_data_frame[i,] <- c(simulated_fractions['lamp5'],simulated_fractions['pvalb'],
                                     simulated_fractions['sst'],simulated_fractions['vip'])
    }
    return(simulated_data_frame)
}

# gets the fraction values
get_pc_fractions <- function(input_df,condition){
    subset_df <- input_df[grep(condition, rownames(input_df)), ]
    subset_df <- subset_df[,1:4]
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
    cluster_list <- c('lamp5','pvalb','sst','vip')
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
    names(z_df)[1] <- paste('lamp5')
    names(z_df)[2] <- paste('pvalb')
    names(z_df)[3] <- paste('sst')
    names(z_df)[4] <- paste('vip')
    return(z_df)
}

mm_z_df <- get_z_df('mm')
run_z_df <- get_z_df('run')
vis_z_df <- get_z_df('vis')
                              
lamp5_df <- cbind(mm_z_df[,1,drop = FALSE],run_z_df[,1,drop = FALSE],vis_z_df[,1,drop = FALSE]) 
colnames(lamp5_df) <- c("Mismatch","Running","Visual")
                              
pvalb_df <- cbind(mm_z_df[,2,drop = FALSE],run_z_df[,2,drop = FALSE],vis_z_df[,2,drop = FALSE]) 
colnames(pvalb_df) <- c("Mismatch","Running","Visual")
                              
sst_df <- cbind(mm_z_df[,3,drop = FALSE],run_z_df[,3,drop = FALSE],vis_z_df[,3,drop = FALSE]) 
colnames(sst_df) <- c("Mismatch","Running","Visual")
                              
vip_df <- cbind(mm_z_df[,4,drop = FALSE],run_z_df[,4,drop = FALSE],vis_z_df[,4,drop = FALSE]) 
colnames(vip_df) <- c("Mismatch","Running","Visual")

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
geom_hline(yintercept=0, col = "black", size = 1) + ylim(-5, 5) +
geom_hline(yintercept=3.196720, col = "gray", size = 1, linetype = 'dotted') +   #bonferoni corrected from 0.05
geom_hline(yintercept=-3.196720, col = "gray", size = 1,linetype = 'dotted') +
theme(axis.text.x = element_text(angle = 70, vjust = 0.97, hjust=1.0)) +
theme(axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm")) +
theme(axis.title.x=element_blank())  +
theme(legend.position="right") +
theme(axis.line = element_line(colour = 'black', size = 1)) +
theme(axis.ticks = element_line(colour = "black", size = 1)) +
theme(legend.position = "none")
return(monte_plot)
}

lamp5_monte <- shiny_monte_plot(c((lamp5_df[,1]),(lamp5_df[,2]),(lamp5_df[,3])))+ ggtitle('lamp5') +
theme(text=element_text(size=30)) + theme(plot.title = element_text(size=30)) +
theme(plot.title = element_text(hjust = 0.5)) +
theme(text=element_text(family="ArialMT"))
                              
pvalb_monte <- shiny_monte_plot(c((pvalb_df[,1]),(pvalb_df[,2]),(pvalb_df[,3])))+ ggtitle('pvalb') +
theme(text=element_text(size=30)) + theme(plot.title = element_text(size=30)) +
theme(plot.title = element_text(hjust = 0.5)) +
theme(text=element_text(family="ArialMT"))
            
sst_monte <- shiny_monte_plot(c((sst_df[,1]),(sst_df[,2]),(sst_df[,3])))+ ggtitle('sst') +
theme(text=element_text(size=30)) + theme(plot.title = element_text(size=30)) +
theme(plot.title = element_text(hjust = 0.5)) +
theme(text=element_text(family="ArialMT"))
                
vip_monte <- shiny_monte_plot(c((vip_df[,1]),(vip_df[,2]),(vip_df[,3])))+ ggtitle('vip') +
theme(text=element_text(size=30)) + theme(plot.title = element_text(size=30)) +
theme(plot.title = element_text(hjust = 0.5)) +
theme(text=element_text(family="ArialMT"))

panel.a <- lamp5_monte
panel.b <- pvalb_monte
panel.c <- sst_monte
panel.d <- vip_monte

panels = list(panel.a,panel.b,panel.c,panel.d)
             
lay <- rbind(c(1,1,2,2,3,3,4,4),
             c(1,1,2,2,3,3,4,4))

fig_path <- get_path(variables$figures,'fig_s3_b_through_e.pdf')
pdf(fig_path, width = 30, height = 10)
grid.arrange(grobs = panels, layout_matrix = lay)
dev.off()