# Single-cell RNA-seqencing analysis of functionally labelled neuronal populations:

This repository contains the single-cell RNA-sequencing code I have written for my publication entitled: [Molecularly targetable cell types in mouse visual cortex have distinguishable prediction error responses](https://www.cell.com/neuron/pdf/S0896-6273(23)00626-8.pdf).

<sub>**Please note** that this repository, even with the appropriate libraries and packages installed, will not operate independently. Due to size limitations, the **original datasets** are not included. However, for those who are interested, the *original published dataset* as well as the code are available [elsewhere](https://doi.org/10.5281/zenodo.8229544).</sub>
<br>

<sub>**Disclaimer**: this README is still under construction.<sub>

## Project Organization
```

┌── R_figs_master.r/                                : pipeline control
├── images/                                         : contains example images used for explanations within the README
│   └── XXX.png
│   └── XXX.png
|── generate_allen_data.r                           : formats refernce dataset
├── preprocess.r/                                   : processes countfiles from cellrnager
├── merge_and_split.r                               : switching sample formats, normalization and merging files into a single object
├──  feature_selection_and_correction.r             : performing feature selection, adding chemistry information,feature exclusion and GLM correction
├── perform_liger.r                                 : non-negative matrix factorization with liger for sample integration
├── calculate_nearest_neighbors.r                   : cell group identification with weighted nearest neigbor
├── filter_sce.r                                    : performing cluster exclusion and distant neighbor exclusion 
├── l23_scmap_assignment.r                          : subtype identification with scmap
├── get_the_figs.r                                  : summary statistics, bootstrap analysis, UMAP plots and figure generation
├── fig_3c.r                                        : correlation heatmap and bootstrap testing for cell type assignment
├── fig_s2_c_through_h.r                            : additional UMAP visualizations
├── fig_s2b.r                                       : UMAP visulisation for reference data set.
├── fig_s3_b_through_e.r                            : bootstrap analysis for enrichment of interneuron groups across photoconversion conditions.
├── campari2_genome_construction.py                 : constructs modified version of mouse genome to account for viral expression
├── single_cell_mapping_pipeline.py                 : used to pass multiple files through the cellranger pipeline
├── single_cell_variables.r                         : stores variables parameters passed to the rest of the pipeline.
├── get_path.r                                      : simple file path construction function
├── LICENSE.md                                      : license
└── README.md                                       : project description

```
<br>



## Genome Construction

```{python}
# adding a custom chromosome
campari2 = ['V','Campari2','exon','1','3500','.','+','.',  ###first 8 elements, tab delimited
            'gene_id "Campari2"',   ### remaining elements, semicolon delimited
            ' gene_version "1"',
            ' transcript_id "Campari2"',
            ' transcript_version "1"',
            ' gene_name "Campari2"',
            ' gene_source "addgene"',
            ' gene_biotype "protein_coding"',
            ' transcript_name "Campari2"',
            ' transcript_source "custom"',
            ' transcript_biotype "protien_coding"',
            ' protein_id "Campari2"',
            ' tag "basic"']
```

Initial processing of the single-cell RNA-sequencing data, up until the construction of the cell-barcode count matrices, of the single-cell RNA-sequencing data was done with CellRanger (3.1.0). Reads were mapped against a custom genome constructed from the Genome Reference Consortium Mouse Build 39 with a version 104 GTF file. Viral expression was accounted for by appending the CaMPARI2 plasmid sequence (Addgene) to the end of the FASTA file and the GTF file was edited to include the regions present on the viral mRNA. GTF and FASTA edits were performed in Python, and the final genome build was constructed with CellRanger. All subsequent analyses on the raw cell feature barcode matrices were done in R. Expression data was imported with the SCATER package.

## Mapping and preprocessing the data

<p align="center">
            <img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/droplet_utils_example.png" height= "600" align = "middle">
</p>

Expression data was imported with the SCATER package. 59 Cells were distinguished from empty droplets with DropletUTILS, 60 using an FDR of 0.05% and a lower limit for cell detection that was calculated on a per dataset basis or set to the median UMI count for droplets greater than 10 UMIs. Cells with a mitochondrial read percentage greater than 20% as well as cells with read counts that were more than three standard deviations below the sample average were excluded. Feature selection and UMI correction were done with SCtransform from the Seurat package. 61,62 Initially, 2000 genes were selected and mitochondrial percentage, sample effects, 10x chemistry (v2 vs. v3), and total feature counts were regressed out. Selected genes were filtered to remove any genes that were differentially expressed between groups of cells separated by read count or mitochondrial percentage, thus reducing the impact of cell quality on the analysis. The final corrected UMI matrix was then used for alignment with a reference dataset. 

##  Reference data set preparation

```{r}
generate_allen_data <- function(){
allen_path = 'data_repository/single_cell_data/allen_data/mouse_VISp_2018-06-14_exon-matrix.csv'
allen.rna <- as.sparse(read.csv(file = allen_path, sep = ",", header = TRUE, row.names = 1))
gene_path <- 'data_repository/single_cell_data/allen_data/mouse_VISp_2018-06-14_genes-rows.csv'
gene.names <- read.csv(file = gene_path, header = TRUE, sep = ',' )
gene_symbol <- gene.names[ ,'gene_symbol']
rownames(allen.rna) <- gene_symbol
meta_path <-  'data_repository/single_cell_data/allen_data/mouse_VISp_2018-06-14_samples-columns.csv'
sample.info <- read.csv(file = meta_path, header = TRUE, sep = ',' )
allen_cluster <- sample.info[ ,'cluster']
allen_cluster <- sapply(allen_cluster, toString)
allen_sub_class <- sample.info[ ,'subclass']
allen_sub_class <- sapply(allen_sub_class, toString)
allen_class <- sample.info[ ,'class']
allen_class <- sapply(allen_class, toString)
allen_region <- sample.info[ ,'brain_region']
allen_region <- sapply(allen_region, toString)
allen_subregion <- sample.info[ ,'brain_subregion']
allen_subregion <- sapply(allen_subregion, toString)
allen.cells <- CreateSeuratObject(counts = allen.rna)
allen.cells@meta.data$allenCluster <- allen_cluster
allen.cells@meta.data$allenSubClass <- allen_sub_class
allen.cells@meta.data$allenClass <- allen_class
allen.cells@meta.data$allenRegion <- allen_region
allen.cells@meta.data$allenSubregion <- allen_subregion
allen.cells@meta.data$sample <- 'allen'
allen.cells@meta.data$facs_sample <- 'allen'
allen.cells@meta.data$Batch <- 'allen'
allen.cells@meta.data$technology <- 'smartseq'
allen.cells@meta.data$percent.mt <- sample.info$percent_mt_exon_reads
allen.cells@meta.data$percent.rb <- sample.info$percent_rrna_reads
allen.cells@meta.data$percent.virus <- 0
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
layers.list <- c("L1","L2/3","L4","L1-L4","L1-L2/3","L2/3-L4")
cluster.list <- c('L2/3 IT VISp Adamts2','L2/3 IT VISp Rrad','L2/3 IT VISp Agmat','L4 IT VISp Rspo1')
exclusion.list <- c('High Intronic','No Class', 'Low Quality','Batch Grouping','Doublet','L5 IT','L6 IT','L5 PT','Endo','CR','NP')
allen.cells <- subset(x = allen.cells, subset = allenCluster %in% cluster.list 
               | allenSubregion %in% layers.list
               & allenSubClass %not in% exclusion.list)
allen.cells <- NormalizeData(allen.cells)
allen.cells <- FindVariableFeatures(allen.cells, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
allen.cells <- ScaleData(allen.cells)
allen.cells <- RunPCA(allen.cells, verbose = FALSE)
allen.cells <- RunUMAP(allen.cells, reduction = "pca", dims = 1:25, verbose = FALSE)
Idents(object = allen.cells) <- 'allenSubClass'
saveRDS(allen.cells, 'stored_variables/allen.cells.rds')
```

To align cells from this study with the reference dataset some portions of the data were excluded that would not have been relevant for the dataset I collected.....

## Data set integration and type assignment 

To align cells from this study with the reference dataset, 25 excluding L5 and L6 neurons, we used the Liger package 34 with a k-value of 12 and lambda-value of 10. These parameters were chosen based on how well groups were separated on a UMAP plot, as well as the expression of various markers normally enriched in each major cortical cell group (Figure S2A). Then, iNMF factors were used to assign each cell from this study to a cell group defined in the reference dataset (L4, L2/3, Sst etc.) using a weighted nearest neighbor algorithm. For each cell in our dataset, we determined the 10 nearest neighbors in the reference dataset, and each of these 10 neighbors were assigned a value corresponding to the inverse of the distance to the cell in our dataset. The 10 nearest neighbors were pooled by their cell group identity, and the inverse distance values summed across cell groups. Cell identity of our cell was defined as the cell group with the highest summed inverse distance value. The initial cell total was 65 002. However, we applied several filters to the dataset. First, since the initially assigned Pvalb group contained many cells that expressed excitatory neuron markers graphical clustering was used to subselect this group and isolate Pvalb positive neurons that were negative for excitatory markers. This eliminated 4762 cells from the dataset. Second, cells whose distance to the closest cell in the reference dataset exceeded 0.002 in the iNMF factor space were removed from the dataset, a filter that eliminated 9724 cells. Third, we excluded cell types that con-
tained less than 300 cells across the entire dataset (eliminating 67 cells). In total 2353 (10.07%) L2/3 neurons were excluded from further analysis. The final cell count was 50 449.  

<p align="center">
<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/dotplot_supp.png" height= 600>
</p>

For subtype assignment within L2/3: Adamts2, Rrad, and Agmat cell type assign-ments were determined separately within the L2/3 group with scMAP 63 using the scmapCluster function; using scMAP’s selectFeatures function, 250 features (genes) were initially selected from the reference dataset, 8 of which were not present in our 10x dataset. This procedure uniquely assigned a cell type identity to all but 60 of the cells that remained unassigned and were excluded from further analysis. Unassigned cells had differing maximum cosine similarities, maximum correlations, or maximum rank order correlations with different L2/3 excitatory cell types. Thus, 242 genes were used (with no similarity cutoff) for cell type assignment. These parameters allowed for the assignment of L2/3 neurons to their corresponding cell types in the reference dataset 25 types (Figure 3). Expression values used in Figure 3 were adjusted with the ScaleData function from the Seurat package and regressed against: UMI count, mitochondrial percentage, sample identity, sequencing chemistry, and total feature count to remove unwanted sources of variance.

<p align="center">
<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/heatmaps_main.png" height=600>
</p>

To quantify the gene expression similarity in the three cell types between our dataset and the reference dataset (Figure 3C), we used a bootstrap estimate of the similarity expected by chance and compared the actual similarity to this bootstrap estimate. We first computed the correlation between the gene expression vectors of all cells of a given cell type with the corresponding average expression vector of the same cell type in the reference dataset. We constructed the bootstrap estimate of random correlation by shuffling cell type identity and repeating the correlation analysis.

<p align="center">
<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/umap_main.png" height= 700>
</p>

## Bootstrap analysis of cell type enrichment

 To analyze how the proportion of L2/3 cell types changed as a func-
tion of photoconversion (Figures 4B–4E), we generated a bootstrap distribution of the fraction of cells we would expect to be of a
given type for all three cell types. We generated 100 000 bootstrap samples of size N, where N was the average number of L2/3 neu-
rons per single-cell sample. This bootstrap distribution was based on the pooled data of all photoconversion groups. We then
compared the actual fraction of each cell type observed in the data for each photoconversion group to the one we would expect
to find by chance. We quantified the difference between the chance distribution and the actual value as a Z score.

<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/bootstrap_main.png" height=400>
