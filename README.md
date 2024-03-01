# Single-cell RNA-seqencing analysis of functionally labelled neuronal populations:

Single-cell RNA-sequencing analysis for the publication entitled: [Molecularly targetable cell types in mouse visual cortex have distinguishable prediction error responses](https://www.cell.com/neuron/pdf/S0896-6273(23)00626-8.pdf).



**Please note** that this repository, even with the appropriate libraries and packages installed, will not operate independently. Due to size limitations, the **original datasets** are not included. However, for those who are interested, the *original published dataset* as well as the code are available [elsewhere](https://doi.org/10.5281/zenodo.8229544).

At the moment this README is still under construction, more details to follow.


## Project Organization
```
.
├── R_figs_master.r/                                : pipeline control
├── Images/                                         : contains example images used for explanations within the README
│   └── XXX.png
│   └── XXX.png
|── generate_allen_data.r                           :formats refernce dataset
├── preprocess.r/                                   :processes countfiles from cellrnager
├── merge_and_split.r    : switching sample formats, normalization and merging files into a single object
├──  feature_selection_and_correction.r             : Performing feature selection, adding chemistry information,feature exclusion and GLM correction
├── perform_liger.r                                 : Non-negative matrix factorization with liger for sample integration
├── calculate_nearest_neighbors.r                   : Cell group identification with weighted nearest neigbor
├── filter_sce.r                                    : Performing cluster exclusion and distant neighbor exclusion 
├── perform_scmap.r                                 : subtype identification with scmap
├── get_the_figs.r                                  : Summary statistics, bootstrap analysis, UMAP plots and figure generation
├── fig_3c.f                                        : XXX
├── fig_s2_c_through_h.r                            : XXX
├── fig_s2b.r                                       : XXX
├── sig_s3_b_through_e.r                            : XX
├── campari2_genome_construction.py                 : Constructs modified version of mouse genome to account for viral expression
├── single_cell_mapping_pipeline.py                 : Used to pass multiple files through the cellranger pipeline
├── LICENSE.md                                      : MIT License
└── README.md                                       : Project description


