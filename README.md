# Single-cell RNA-seqencing analysis of functionally labelled neuronal populations:

This repository contains the single-cell RNA-sequencing code I have written for my publication entitled: [Molecularly targetable cell types in mouse visual cortex have distinguishable prediction error responses](https://www.cell.com/neuron/pdf/S0896-6273(23)00626-8.pdf). This README contains modified excerpts from that study as well as descriptions for each individual piece of code.

<sub>**Please note** that this repository, even with the appropriate libraries and packages installed, will not operate independently. Due to size limitations, the **original datasets** are not included. However, for those who are interested, the *original published dataset* as well as the code are available [elsewhere](https://doi.org/10.5281/zenodo.8229544).</sub>
<br>

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

## Modified excerpts (figures and methods) from O'Toole et al. 2023 relevant to this repository

<p align="center">
<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/umap_main.png" height= 700>
</p>

### Single-cell RNA-sequencing of photoconverted cells
**(A)** Uniform manifold approximation and projection (UMAP) of neurons isolated from the superficial cortical layers (L1, L2/3, and L4), based on the integrative non-negative matrix factorization (iNMF) vectors calculated from a shared projection of our data and a reference dataset.25 Major cell groups are labeled, and inferred cell type identity of L2/3 cells are color coded. Plot order is randomized, except for the L2/3 Rrad cell type that is plotted in the foreground due to its relative scarcity.
**(B)** The L2/3 cluster with cells expressing one or more UMIs corresponding to Adamts2 highlighted. Note, the expression of Adamts2 overlaps with the Adamts2 cell type shown in (A), but the Adamts2 cell type is defined by the expression of a set of different genes of which Adamts2 is just one.
**(C)** As in **(B)**, but for Agmat.
**(D)** As in **(B)**, but for Rrad.

***

<p align="center">
<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/heatmaps_main.png" height=600>
</p>

### L2/3 excitatory cell type assignment
**(A)** Expression heatmap for a select number of genes enriched in each of the three L2/3 excitatory cell types for all L2/3 excitatory neurons examined in this study as in Figure 2A. Unassigned neurons (60) are not shown. Color values correspond to expression level. To reduce aliasing, data are smoothed across cells with a window size of 20.
**(B)** As in **(A)**, but for all L2/3 V1 excitatory neurons from a reference dataset.
**(C)** Marker gene expression of L2/3 excitatory neurons (19,689) examined in this study (6,495 Adamts2, 12,143 Agmat, 991 Rrad) correlated with marker gene expression of all L2/3 V1 excitatory neurons from a reference dataset.25 Correlation of the gene expression vectors was calculated between individual neurons of this study with averaged expression vectors for each of the L2/3 excitatory types of the reference dataset (see STAR Methods). Assigned cell types of each neuron are indicated with color bars on top (this study) and left.25 To prevent aliasing, data are smoothed across cells with a window size of 20.

***

<p align="center">
<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/bootstrap_main.png" height=400>
</p>

### Differential enrichment of L2/3 excitatory transcriptional cell types in different photoconversion groups
**(A)** Fold change in the proportion of L2/3 excitatory cell types between the high and low photoconversion for photoconversion on mismatch, running onset, and grating onset. The total percentages of the three cell types were: Adamts2 33%, Agmat 62%, and Rrad 5.0%. For mismatch photoconversion the percentages in the low and high photoconversion groups were: Adamts2 low 22% and high 36%, Agmat low 73% and high 60%, and Rrad low = 4.9% and high 4.5%. For running onset photoconversion, the percentages in the low and high photoconversion groups were Adamts2 low 41% and high 36%, Agmat low 53% and high 59%, and Rrad low = 6.1% and high 5.3%. For grating onset photoconversion, the percentages in the low and high photoconversion groups were Adamts2 low 41% and high 28%, Agmat low 55% and high 65%, and Rrad low = 3.8% and high 6.9%.
**(B)** Enrichment of the L2/3 Adamts2 cell type in the high (red), intermediate (pink), and low (green) photoconversion groups for the visuomotor mismatch photoconversion experiments, relative to the distribution observed in random samples (Z score distribution shown in black, dashed gray lines correspond to p = 0.01). The L2/3 Adamts2 cell type is significantly enriched in the high visuomotor mismatch photoconversion group.
**(C)** Enrichment for the L2/3 Adamts2 cell type across photoconversion and labeling conditions. Z score values are displayed for visuomotor mismatch (mismatch; as in A), running onset (Running), and grating onset when stationary (Grating) photoconversion experiments.
**(D)** As in **(C)**, but for the L2/3 Agmat cell type.
**(E)** As in **(C)**, but for the L2/3 Rrad cell type.

***

<p align="center">
<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/figure_s2.png" height= 600>
</p>

### Cellular responses to mismatch, running and grating onsets compared to their fluorescence ratio after photoconversion during mismatch.

**(A)** Dot plot of selected marker genes for each of the major cell groups in Figure 2A (with a cumulative population size greater than 300 cells across all samples). Dot size represents the percentage of cells expressing each marker while the color represents the average normalized expression value. Exc.: excitatory neurons, Inh.: inhibitory neurons, Oligo.: Oligodendrocytes, IT: intratelencephalic neurons. 
**(B)** UMAP plot of L2/3 neurons in the reference dataset (Tasic et al., 2018) with the three L2/3 cell types (Adamts2, Agmat, and Rrad) marked as they were defined there.  
**(C)** A UMAP plot of L2/3 excitatory neurons collected using photoconversion on mismatch. Low, intermediate, and high photoconversion groups are labelled as green, gray and red, respectively.  
**(D)** As in **B**, but for photoconversion during gratings. 
**(E)** A density plot for the high photoconversion group collected using photoconversion on mismatch, calculated based on the UMAP plot shown in **B**. 
**(F)** A density plot for the high photoconversion group collected using photoconversion on gratings, calculated based on the UMAP plot shown in **C**. 
**(G)** A density plot of cells expressing at least one copy of Adamts2 (data as in Figure 2B). 
**(H)** As in **E**, but for Baz1a. 

***

<p align="center">
<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/figure_s3.png" height= 600>
</p>

### Abundance of major cell groups and L2/3 excitatory types and enrichment analysis for inhibitory neuron types. Related to Figures 2-4. 

**(A)** Abundance of each major cell group as shown in Figure 2A, across the three photoconversion groups and the three photoconversion types. mm: Mismatch; run: Running onsets; vis: Grating onsets; 1: Low photoconversion group; 2: Intermediate photoconversion group; 3: High photoconversion group. 
**(B)** Enrichment for the Lamp5 cell type across photoconversion and labeling conditions. Z-score values are displayed for visuomotor mismatch (mismatch), running onset (running), and grating onset when stationary (grating) labeling experiments. Enrichment is calculated relative to the major inhibitory neuron groups only. 
**(C)** As in **B**, but for the Pvalb cell type. 
**(D)** As in **B**, but for the Sst cell type. 
**(E)** As in **B**, but for the Vip cell type.  

***

### Methods

Initial processing of the single-cell RNA-sequencing data, up until the construction of the cell-barcode count matrices, of the single-cell RNA-sequencing data was done with CellRanger (3.1.0). Reads were mapped against a custom genome constructed from the Genome Reference Consortium Mouse Build 39 with a version 104 GTF file. Viral expression was accounted for by appending the CaMPARI2 plasmid sequence (Addgene) to the end of the FASTA file and the GTF file was edited to include the regions present on the viral mRNA. GTF and FASTA edits were performed in Python, and the final genome build was constructed with CellRanger. All subsequent analyses on the raw cell feature barcode matrices were done in R.

Expression data was imported with the SCATER package. 59 Cells were distinguished from empty droplets with DropletUTILS,60 using an FDR of 0.05% and a lower limit for cell detection that was calculated on a per dataset basis or set to the median UMI count for droplets greater than 10 UMIs. Cells with a mitochondrial read percentage greater than 20% as well as cells with read counts that were more than three standard deviations below the sample average were excluded. Feature selection and UMI correction were done with SCtransform from the Seurat package.61,62 Initially, 2000 genes were selected and mitochondrial percentage, sample effects, 10x chemistry (v2 vs. v3), and total feature counts were regressed out. Selected genes were filtered to remove any genes that were differentially expressed between groups of cells separated by read count or mitochondrial percentage, thus reducing the impact of cell quality on the analysis. The final corrected UMI matrix was then used for alignment with a reference dataset.

To align cells from this study with the reference dataset,25 excluding L5 and L6 neurons, we used the Liger package34 with a k-value of 12 and lambda-value of 10. These parameters were chosen based on how well groups were separated on a UMAP plot, as well as the expression of various markers normally enriched in each major cortical cell group. Then, iNMF factors were used to assign each cell from this study to a cell group defined in the reference dataset (L4, L2/3, Sst etc.) using a weighted nearest neighbor algorithm. For each cell in our dataset, we determined the 10 nearest neighbors in the reference dataset, and each of these 10 neighbors were assigned a value corresponding to the inverse of the distance to the cell in our dataset. The 10 nearest neighbors were pooled by their cell group identity, and the inverse distance values summed across cell groups. Cell identity of our cell was defined as the cell group with the highest summed inverse distance value. The initial cell total was 65 002. However, we applied several filters to the dataset. First, since the initially assigned Pvalb group contained many cells that expressed excitatory neuron markers, graphical clustering was used to subselect this group and isolate Pvalb positive neurons that were negative for excitatory markers. This eliminated 4762 cells from the dataset. Second, cells whose distance to the closest cell in the reference dataset exceeded 0.002 in the iNMF factor space were removed from the dataset, a filter that eliminated 9724 cells. Third, we excluded cell types that contained less than 300 cells across the entire dataset (eliminating 67 cells). In total 2353 (10.07%) L2/3 neurons were excluded from further analysis. The final cell count was 50 449. For subtype assignment within L2/3: Adamts2, Rrad, and Agmat cell type assignments were determined separately within the L2/3 group with scMAP63 using the scmapCluster function; using scMAP’s selectFeatures function, 250 features (genes) were initially selected from the reference dataset, 8 of which were not present in our 10x dataset. This procedure uniquely assigned a cell type identity to all but 60 of the cells that remained unassigned and were excluded from further analysis. Unassigned cells had differing maximum cosine similarities, maximum correlations, or maximum rank order correlations with different L2/3 excitatory cell types. Thus, 242 genes were used (with no similarity cutoff) for cell type assignment. These parameters allowed for the assignment of L2/3 neurons to their corresponding cell types in the reference dataset25 types. Expression values used in the heatmap figure were adjusted with the ScaleData function from the Seurat package and regressed against: UMI count, mitochondrial percentage, sample identity, sequencing chemistry, and total feature count to remove unwanted sources of variance. To quantify the gene expression similarity in the three cell types between our dataset and the reference dataset, we used a bootstrap estimate of the similarity expected by chance and compared the actual similarity to this bootstrap estimate. We first computed the correlation between the gene expression vectors of all cells of a given cell type with the corresponding average expression vector of the same cell type in the reference dataset. We constructed the bootstrap estimate of random correlation by shuffling cell type identity and repeating the correlation analysis. To analyze how the proportion of L2/3 cell types changed as a function of photoconversion, we generated a bootstrap distribution of the fraction of cells we would expect to be of a given type for all three cell types. We generated 100 000 bootstrap samples of size N, where N was the average number of L2/3 neurons per single-cell sample. This bootstrap distribution was based on the pooled data of all photoconversion groups. We then compared the actual fraction of each cell type observed in the data for each photoconversion group to the one we would expect to find by chance. We quantified the difference between the chance distribution and the actual value as a Z score.


