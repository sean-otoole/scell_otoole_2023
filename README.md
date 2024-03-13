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

## Mapping and preprocessing the data

<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/droplet_utils_example.png" height= 600>

##  Reference data set preparation

## Data set integration and visualization 

<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/umap_main.png" height= 700>

## Assessing group assignment

<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/dotplot_supp.png" height= 600>

## Assessing type assignment

<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/heatmaps_main.png" height=600>

## Bootstrap analysis of cell type enrichment

<img src="https://github.com/sean-otoole/scell_otoole_2023/blob/main/images/bootstrap_main.png" height=400>
