# Analysis code overview
### To complete a reproducible run from raw data:
Retrieve the fastq files from SRA and then align to the equine genome (EquCab3.0; Ensembl) using Cell Ranger (version 6.1.2) with default settings. 
Instructions to download the fastq files can be found in [:file\_folder: input](/input). 
The alignment code is currently not provided, but can be shared if desired. If you want it, you can request through GitHub.

### To complete a reproducible run from cellranger output count matricies:
Download the supplementary zip folder on the NCBI GEO project page as decribed in [:file\_folder: input](/input) directory.
In addition to the count matrices you will also need the .csv files located in the [:file\_folder: sheets](/input/sheets) directory.
Once all files are in place, you can use the provided scripts (generally split by major analysis step) and the source file [customFunctions.R](/analysis/customFunctions.R) to explore the analysis approach.

NOTE: the file paths in the indivdual scripts will may need to be modified when reading/writing data to run on your system.

### To reproduce/explore data using processed data:
Download the .rds files on the Zenodo project page as described in [:file\_folder: input](/input) directory.
These files can be loaded into individual analysis scripts to bypass some of the more computationally intense steps, or they can be used for your own exploration!

### Envrionement to reproduce:
This work was completed on the UC Boulder Alpine supercomputer on a Linux operating system. All code was run in a Conda environment and I am happy to share the .yml file upon request.

Otherwise, see below for packages loaded into the system.
### Software versions
```r
> sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.4 (Ootpa)

Matrix products: default
BLAS/LAPACK: /projects/dyammons@colostate.edu/software/anaconda/envs/r_env/lib/libopenblasp-r0.3.18.so

locale:
 [1] LC_CTYPE=C.UTF-8    LC_NUMERIC=C        LC_TIME=C          
 [4] LC_COLLATE=C        LC_MONETARY=C       LC_MESSAGES=C      
 [7] LC_PAPER=C          LC_NAME=C           LC_ADDRESS=C       
[10] LC_TELEPHONE=C      LC_MEASUREMENT=C    LC_IDENTIFICATION=C

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] CellChat_1.5.0              igraph_1.5.1               
 [3] UpSetR_1.4.0                scProportionTest_0.0.0.9000
 [5] ComplexHeatmap_2.13.2       ggtree_3.2.1               
 [7] ape_5.7-1                   scuttle_1.4.0              
 [9] scRNAseq_2.8.0              ggpubr_0.4.0               
[11] slingshot_2.7.0             TrajectoryUtils_1.2.0      
[13] SingleCellExperiment_1.16.0 princurve_2.1.6            
[15] clusterProfiler_4.2.2       msigdbr_7.5.1              
[17] ggsankey_0.0.99999          lemon_0.4.5                
[19] reshape_0.8.9               viridis_0.6.2              
[21] viridisLite_0.4.1           SingleR_1.8.1              
[23] SeuratDisk_0.0.0.9019       RColorBrewer_1.1-3         
[25] pheatmap_1.0.12             DESeq2_1.34.0              
[27] SummarizedExperiment_1.24.0 Biobase_2.54.0             
[29] MatrixGenerics_1.6.0        matrixStats_1.0.0          
[31] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
[33] IRanges_2.28.0              S4Vectors_0.32.4           
[35] BiocGenerics_0.40.0         colorspace_2.0-3           
[37] ggrepel_0.9.1               cowplot_1.1.1              
[39] scales_1.2.1                patchwork_1.1.2            
[41] DoubletFinder_2.0.3         clustree_0.4.4             
[43] ggraph_2.0.5                forcats_0.5.2              
[45] stringr_1.4.1               dplyr_1.0.10               
[47] purrr_0.3.5                 readr_2.1.2                
[49] tidyr_1.2.1                 tibble_3.1.8               
[51] ggplot2_3.3.6               tidyverse_1.3.1            
[53] SeuratObject_4.1.3          Seurat_4.3.0               

loaded via a namespace (and not attached):
  [1] statnet.common_4.7.0          rsvd_1.0.5                   
  [3] svglite_2.1.0                 ica_1.0-3                    
  [5] Rsamtools_2.10.0              foreach_1.5.2                
  [7] lmtest_0.9-40                 crayon_1.5.2                 
  [9] MASS_7.3-58.1                 nlme_3.1-160                 
 [11] backports_1.4.1               reprex_2.0.1                 
 [13] GOSemSim_2.20.0               rlang_1.0.6                  
 [15] XVector_0.34.0                ROCR_1.0-11                  
 [17] readxl_1.4.1                  irlba_2.3.5                  
 [19] filelock_1.0.2                BiocParallel_1.37.1          
 [21] rjson_0.2.21                  bit64_4.0.5                  
 [23] glue_1.6.2                    rngtools_1.5.2               
 [25] sctransform_0.3.5             parallel_4.1.1               
 [27] spatstat.sparse_3.0-0         AnnotationDbi_1.56.2         
 [29] DOSE_3.20.1                   spatstat.geom_3.0-5          
 [31] haven_2.4.3                   tidyselect_1.2.0             
 [33] fitdistrplus_1.1-8            XML_3.99-0.12                
 [35] zoo_1.8-9                     GenomicAlignments_1.30.0     
 [37] xtable_1.8-4                  magrittr_2.0.3               
 [39] cli_3.6.0                     zlibbioc_1.40.0              
 [41] rstudioapi_0.14               miniUI_0.1.1.1               
 [43] sp_1.5-1                      fastmatch_1.1-3              
 [45] ensembldb_2.18.3              treeio_1.18.1                
 [47] shiny_1.7.1                   BiocSingular_1.10.0          
 [49] xfun_0.39                     clue_0.3-62                  
 [51] cluster_2.1.4                 tidygraph_1.2.2              
 [53] KEGGREST_1.34.0               interactiveDisplayBase_1.32.0
 [55] listenv_0.8.0                 Biostrings_2.62.0            
 [57] png_0.1-7                     future_1.29.0                
 [59] withr_2.5.0                   bitops_1.0-7                 
 [61] ggforce_0.4.1                 plyr_1.8.7                   
 [63] cellranger_1.1.0              AnnotationFilter_1.18.0      
 [65] coda_0.19-4                   pillar_1.8.1                 
 [67] GlobalOptions_0.1.2           cachem_1.0.6                 
 [69] GenomicFeatures_1.46.5        fs_1.5.2                     
 [71] hdf5r_1.3.8                   GetoptLong_1.0.5             
 [73] DelayedMatrixStats_1.16.0     vctrs_0.5.1                  
 [75] ellipsis_0.3.2                generics_0.1.3               
 [77] NMF_0.24.0                    tools_4.1.1                  
 [79] munsell_0.5.0                 tweenr_2.0.2                 
 [81] fgsea_1.20.0                  DelayedArray_0.20.0          
 [83] fastmap_1.1.0                 compiler_4.1.1               
 [85] abind_1.4-5                   httpuv_1.6.5                 
 [87] rtracklayer_1.54.0            ExperimentHub_2.2.1          
 [89] pkgmaker_0.32.2.900           plotly_4.10.1                
 [91] GenomeInfoDbData_1.2.7        gridExtra_2.3                
 [93] lattice_0.20-45               deldir_1.0-6                 
 [95] utf8_1.2.2                    later_1.3.0                  
 [97] BiocFileCache_2.2.1           jsonlite_1.8.3               
 [99] ScaledMatrix_1.2.0            tidytree_0.3.9               
[101] pbapply_1.5-0                 carData_3.0-5                
[103] sparseMatrixStats_1.6.0       genefilter_1.76.0            
[105] lazyeval_0.2.2                promises_1.2.0.1             
[107] car_3.1-0                     doParallel_1.0.17            
[109] goftest_1.2-3                 sna_2.7                      
[111] spatstat.utils_3.0-1          reticulate_1.34.0            
[113] Rtsne_0.16                    downloader_0.4               
[115] uwot_0.1.14                   survival_3.4-0               
[117] yaml_2.3.6                    systemfonts_1.0.4            
[119] htmltools_0.5.3               memoise_2.0.1                
[121] BiocIO_1.4.0                  locfit_1.5-9.6               
[123] graphlayouts_0.8.3            digest_0.6.30                
[125] assertthat_0.2.1              mime_0.12                    
[127] rappdirs_0.3.3                registry_0.5-1               
[129] RSQLite_2.2.18                yulab.utils_0.0.5            
[131] future.apply_1.10.0           data.table_1.14.4            
[133] blob_1.2.3                    splines_4.1.1                
[135] AnnotationHub_3.2.2           ProtGenerics_1.26.0          
[137] RCurl_1.98-1.12               broom_1.0.1                  
[139] hms_1.1.2                     modelr_0.1.8                 
[141] BiocManager_1.30.19           shape_1.4.6                  
[143] aplot_0.1.2                   Rcpp_1.0.11                  
[145] RANN_2.6.1                    circlize_0.4.16              
[147] enrichplot_1.14.2             fansi_1.0.3                  
[149] tzdb_0.3.0                    parallelly_1.32.1            
[151] R6_2.5.1                      ggridges_0.5.4               
[153] lifecycle_1.0.3               curl_4.3.3                   
[155] ggsignif_0.6.4                leiden_0.4.3                 
[157] DO.db_2.9                     Matrix_1.5-1                 
[159] qvalue_2.26.0                 RcppAnnoy_0.0.20             
[161] iterators_1.0.14              spatstat.explore_3.0-5       
[163] htmlwidgets_1.5.4             beachmat_2.10.0              
[165] polyclip_1.10-4               biomaRt_2.50.3               
[167] network_1.18.0                shadowtext_0.1.1             
[169] timechange_0.1.1              gridGraphics_0.5-1           
[171] rvest_1.0.3                   globals_0.16.1               
[173] spatstat.random_3.1-3         progressr_0.11.0             
[175] codetools_0.2-18              lubridate_1.9.0              
[177] FNN_1.1.3.1                   GO.db_3.14.0                 
[179] prettyunits_1.1.1             dbplyr_2.1.1                 
[181] RSpectra_0.16-1               gridBase_0.4-7               
[183] gtable_0.3.1                  DBI_1.1.3                    
[185] ggalluvial_0.12.3             ggfun_0.0.8                  
[187] tensor_1.5                    httr_1.4.4                   
[189] KernSmooth_2.23-20            stringi_1.7.8                
[191] progress_1.2.2                reshape2_1.4.4               
[193] farver_2.1.1                  annotate_1.72.0              
[195] xml2_1.3.3                    BiocNeighbors_1.12.0         
[197] restfulr_0.0.15               geneplotter_1.72.0           
[199] ggplotify_0.1.0               scattermore_0.8              
[201] BiocVersion_3.14.0            bit_4.0.4                    
[203] scatterpie_0.1.7              spatstat.data_3.0-0          
[205] pkgconfig_2.0.3               babelgene_22.9               
[207] rstatix_0.7.0                 knitr_1.40  
```
