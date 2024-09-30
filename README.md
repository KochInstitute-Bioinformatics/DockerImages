# Tutorial: how to build a Docker image

Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)

Create an account. Remember your personnal id, this will be important later (mine for example is yannvrb56).

Create a directory, one for each docker image. Avoid spaces, capital letters and special characters for the name of the directory.
For example, let's create a directory for bulk RNA-seq:

```sh
PS C:\Users\yannvrb\Documents\Docker files\bulkrnaseq
```

Inside, create a file named Dockerfile.
Use a text editor to implement commands to install packages.
Save your Dockerfile.

Run docker build command to build your image. Run this command inside bulkrnaseq directory. Replace my personnal id with yours:

```sh
docker build . -t yannvrb56/bulkrnaseq
```

Then push the image to your docker account. Replace my personnal id with yours:

```sh
docker push yannvrb56/bulkrnaseq
```

# Docker Images and their content

## Bulk RNA-seq

To pull the image from the command line:

```sh
docker pull yannvrb56/bulkrnaseq
```

To use as a base image in a new Dockerfile:

```sh
FROM yannvrb56/bulkrnaseq
```

Dockerfile can be found:

```sh
Dockerfiles/bulkrnaseq/Dockerfile
```

Session info:
```sh
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] gplots_3.1.3.1              readxl_1.4.3               
 [3] UpSetR_1.4.0                fgsea_1.30.0               
 [5] cluster_2.1.6               gprofiler2_0.2.3           
 [7] edgeR_4.2.0                 limma_3.60.3               
 [9] tximport_1.32.0             ComplexHeatmap_2.20.0      
[11] apeglm_1.26.1               DESeq2_1.44.0              
[13] SummarizedExperiment_1.34.0 Biobase_2.64.0             
[15] MatrixGenerics_1.16.0       rtracklayer_1.64.0         
[17] GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
[19] IRanges_2.38.1              S4Vectors_0.42.1           
[21] BiocGenerics_0.50.0         ggrepel_0.9.5              
[23] XML_3.99-0.17               matrixStats_1.3.0          
[25] reprex_2.1.0                lubridate_1.9.3            
[27] forcats_1.0.0               stringr_1.5.1              
[29] dplyr_1.1.4                 purrr_1.0.2                
[31] readr_2.1.5                 tidyr_1.3.1                
[33] tibble_3.2.1                ggplot2_3.5.1              
[35] tidyverse_2.0.0             openxlsx_4.2.5.2           

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3       rstudioapi_0.16.0        jsonlite_1.8.8          
  [4] shape_1.4.6.1            magrittr_2.0.3           farver_2.1.2            
  [7] rmarkdown_2.27           GlobalOptions_0.1.2      fs_1.6.4                
 [10] BiocIO_1.14.0            zlibbioc_1.50.0          vctrs_0.6.5             
 [13] Rsamtools_2.20.0         RCurl_1.98-1.14          htmltools_0.5.8.1       
 [16] S4Arrays_1.4.1           curl_5.2.1               cellranger_1.1.0        
 [19] SparseArray_1.4.8        sass_0.4.9               KernSmooth_2.23-24      
 [22] bslib_0.7.0              htmlwidgets_1.6.4        plyr_1.8.9              
 [25] plotly_4.10.4            cachem_1.1.0             GenomicAlignments_1.40.0
 [28] mime_0.12                lifecycle_1.0.4          iterators_1.0.14        
 [31] pkgconfig_2.0.3          Matrix_1.7-0             R6_2.5.1                
 [34] fastmap_1.2.0            shiny_1.8.1.1            GenomeInfoDbData_1.2.12 
 [37] clue_0.3-65              digest_0.6.35            numDeriv_2016.8-1.1     
 [40] colorspace_2.1-0         crosstalk_1.2.1          labeling_0.4.3          
 [43] fansi_1.0.6              timechange_0.3.0         httr_1.4.7              
 [46] abind_1.4-5              compiler_4.4.1           bit64_4.0.5             
 [49] withr_3.0.0              doParallel_1.0.17        BiocParallel_1.38.0     
 [52] highr_0.11               MASS_7.3-60.2            DelayedArray_0.30.1     
 [55] rjson_0.2.21             caTools_1.18.2           gtools_3.9.5            
 [58] tools_4.4.1              httpuv_1.6.15            zip_2.3.1               
 [61] glue_1.7.0               restfulr_0.0.15          promises_1.3.0          
 [64] generics_0.1.3           gtable_0.3.5             tzdb_0.4.0              
 [67] data.table_1.15.4        hms_1.1.3                utf8_1.2.4              
 [70] XVector_0.44.0           foreach_1.5.2            pillar_1.9.0            
 [73] vroom_1.6.5              emdbook_1.3.13           later_1.3.2             
 [76] circlize_0.4.16          lattice_0.22-6           bit_4.0.5               
 [79] tidyselect_1.2.1         locfit_1.5-9.10          Biostrings_2.72.1       
 [82] knitr_1.48               gridExtra_2.3            xfun_0.44               
 [85] statmod_1.5.0            stringi_1.8.4            UCSC.utils_1.0.0        
 [88] lazyeval_0.2.2           yaml_2.3.8               evaluate_0.24.0         
 [91] codetools_0.2-20         bbmle_1.0.25.1           cli_3.6.2               
 [94] xtable_1.8-4             munsell_0.5.1            jquerylib_0.1.4         
 [97] Rcpp_1.0.12              coda_0.19-4.1            png_0.1-8               
[100] bdsmatrix_1.3-7          parallel_4.4.1           bitops_1.0-7            
[103] viridisLite_0.4.2        mvtnorm_1.2-5            scales_1.3.0            
[106] crayon_1.5.2             GetoptLong_1.0.5         rlang_1.1.4             
[109] cowplot_1.1.3            fastmatch_1.1-4         
```

## Single cell RNA-seq and single cell ATAC-seq (Seurat5.1 and Signac)

To pull the image from the command line:

```sh
docker pull yannvrb56/r441seurat5signac
```

To use as a base image in a new Dockerfile:

```sh
FROM yannvrb56/r441seurat5signac
```

Dockerfile can be found:

```sh
Dockerfiles/r441seurat5signac/Dockerfile
```

Dockerfile can be found:

```sh
Dockerfiles/scRNAseq_R422_Seurat43/Dockerfile
```

Session info:
```sh
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.3.0                    TFBSTools_1.42.0                   JASPAR2020_0.99.10                
 [4] BSgenome.Mmusculus.UCSC.mm10_1.4.3 BSgenome_1.72.0                    BiocIO_1.14.0                     
 [7] Biostrings_2.72.1                  XVector_0.44.0                     EnsDb.Mmusculus.v79_2.99.0        
[10] ensembldb_2.28.1                   AnnotationFilter_1.28.0            GenomicFeatures_1.56.0            
[13] AnnotationDbi_1.66.0               Signac_1.14.0                      fgsea_1.30.0                      
[16] cluster_2.1.6                      gprofiler2_0.2.3                   edgeR_4.2.1                       
[19] limma_3.60.4                       ComplexHeatmap_2.20.0              apeglm_1.26.1                     
[22] DESeq2_1.44.0                      rtracklayer_1.64.0                 ggrepel_0.9.6                     
[25] XML_3.99-0.17                      reprex_2.1.0                       lubridate_1.9.3                   
[28] forcats_1.0.0                      stringr_1.5.1                      purrr_1.0.2                       
[31] readr_2.1.5                        tidyr_1.3.1                        tibble_3.2.1                      
[34] tidyverse_2.0.0                    readxl_1.4.3                       ggpubr_0.6.0                      
[37] SingleR_2.6.0                      SummarizedExperiment_1.34.0        Biobase_2.64.0                    
[40] GenomicRanges_1.56.1               GenomeInfoDb_1.40.1                IRanges_2.38.1                    
[43] S4Vectors_0.42.1                   BiocGenerics_0.50.0                MatrixGenerics_1.16.0             
[46] matrixStats_1.4.1                  openxlsx_4.2.7                     ggplot2_3.5.1                     
[49] dplyr_1.1.4                        Seurat_5.1.0                       SeuratObject_5.0.2                
[52] sp_2.1-4                          

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2           poweRlaw_0.80.0             goftest_1.2-3               HDF5Array_1.32.1           
  [5] vctrs_0.6.5                 spatstat.random_3.3-1       digest_0.6.37               png_0.1-8                  
  [9] shape_1.4.6.1               gypsum_1.0.1                deldir_2.0-4                parallelly_1.38.0          
 [13] MASS_7.3-60.2               reshape2_1.4.4              httpuv_1.6.15               foreach_1.5.2              
 [17] withr_3.0.1                 xfun_0.47                   survival_3.6-4              memoise_2.0.1              
 [21] gtools_3.9.5                zoo_1.8-12                  GlobalOptions_0.1.2         pbapply_1.7-2              
 [25] R.oo_1.26.0                 KEGGREST_1.44.1             promises_1.3.0              httr_1.4.7                 
 [29] rstatix_0.7.2               restfulr_0.0.15             globals_0.16.3              fitdistrplus_1.2-1         
 [33] rhdf5filters_1.16.0         rhdf5_2.48.0                rstudioapi_0.16.0           UCSC.utils_1.0.0           
 [37] miniUI_0.1.1.1              generics_0.1.3              curl_5.2.2                  zlibbioc_1.50.0            
 [41] ScaledMatrix_1.12.0         polyclip_1.10-7             GenomeInfoDbData_1.2.12     ExperimentHub_2.12.0       
 [45] SparseArray_1.4.8           pracma_2.4.4                xtable_1.8-4                doParallel_1.0.17          
 [49] evaluate_0.24.0             S4Arrays_1.4.1              BiocFileCache_2.12.0        hms_1.1.3                  
 [53] irlba_2.3.5.1               colorspace_2.1-1            filelock_1.0.3              ROCR_1.0-11                
 [57] reticulate_1.39.0           spatstat.data_3.1-2         magrittr_2.0.3              lmtest_0.9-40              
 [61] later_1.3.2                 lattice_0.22-6              spatstat.geom_3.3-2         future.apply_1.11.2        
 [65] scattermore_1.2             cowplot_1.1.3               RcppAnnoy_0.0.22            pillar_1.9.0               
 [69] nlme_3.1-164                pwalign_1.0.0               iterators_1.0.14            caTools_1.18.3             
 [73] compiler_4.4.1              beachmat_2.20.0             RSpectra_0.16-2             stringi_1.8.4              
 [77] tensor_1.5                  GenomicAlignments_1.40.0    plyr_1.8.9                  crayon_1.5.3               
 [81] abind_1.4-8                 emdbook_1.3.13              locfit_1.5-9.10             bit_4.0.5                  
 [85] fastmatch_1.1-4             codetools_0.2-20            BiocSingular_1.20.0         alabaster.ranges_1.4.2     
 [89] GetoptLong_1.0.5            plotly_4.10.4               mime_0.12                   splines_4.4.1              
 [93] circlize_0.4.16             Rcpp_1.0.13                 fastDummies_1.7.4           dbplyr_2.5.0               
 [97] sparseMatrixStats_1.16.0    cellranger_1.1.0            knitr_1.48                  blob_1.2.4                 
[101] utf8_1.2.4                  seqLogo_1.70.0              clue_0.3-65                 BiocVersion_3.19.1         
[105] fs_1.6.4                    listenv_0.9.1               DelayedMatrixStats_1.26.0   ggsignif_0.6.4             
[109] Matrix_1.7-0                statmod_1.5.0               tzdb_0.4.0                  pkgconfig_2.0.3            
[113] tools_4.4.1                 cachem_1.1.0                RSQLite_2.3.7               viridisLite_0.4.2          
[117] DBI_1.2.3                   numDeriv_2016.8-1.1         celldex_1.14.0              fastmap_1.2.0              
[121] rmarkdown_2.28              scales_1.3.0                ica_1.0-3                   Rsamtools_2.20.0           
[125] broom_1.0.6                 AnnotationHub_3.12.0        coda_0.19-4.1               BiocManager_1.30.25        
[129] dotCall64_1.1-1             carData_3.0-5               RANN_2.6.2                  alabaster.schemas_1.4.0    
[133] farver_2.1.2                yaml_2.3.10                 cli_3.6.3                   leiden_0.4.3.1             
[137] lifecycle_1.0.4             uwot_0.2.2                  mvtnorm_1.3-1               backports_1.5.0            
[141] annotate_1.82.0             BiocParallel_1.38.0         timechange_0.3.0            gtable_0.3.5               
[145] rjson_0.2.23                ggridges_0.5.6              progressr_0.14.0            parallel_4.4.1             
[149] jsonlite_1.8.8              RcppHNSW_0.6.0              bitops_1.0-8                bit64_4.0.5                
[153] Rtsne_0.17                  alabaster.matrix_1.4.2      spatstat.utils_3.1-0        zip_2.3.1                  
[157] CNEr_1.40.0                 bdsmatrix_1.3-7             alabaster.se_1.4.1          R.utils_2.12.3             
[161] spatstat.univar_3.0-1       lazyeval_0.2.2              alabaster.base_1.4.2        shiny_1.9.1                
[165] htmltools_0.5.8.1           GO.db_3.19.1                sctransform_0.4.1           rappdirs_0.3.3             
[169] glue_1.7.0                  TFMPvalue_0.0.9             spam_2.10-0                 httr2_1.0.1                
[173] RCurl_1.98-1.16             gridExtra_2.3               igraph_2.0.3                R6_2.5.1                   
[177] RcppRoll_0.3.1              bbmle_1.0.25.1              Rhdf5lib_1.26.0             DirichletMultinomial_1.46.0
[181] DelayedArray_0.30.1         tidyselect_1.2.1            ProtGenerics_1.36.0         car_3.1-2                  
[185] future_1.34.0               rsvd_1.0.5                  munsell_0.5.1               KernSmooth_2.23-24         
[189] data.table_1.16.0           htmlwidgets_1.6.4           RColorBrewer_1.1-3          rlang_1.1.4                
[193] spatstat.sparse_3.1-0       spatstat.explore_3.3-2      fansi_1.0.6
```

## Single cell RNA-seq (Seurat4.3)

```sh
docker pull yannvrb56/scrnaseq_r422_seurat43
```

To use as a base image in a new Dockerfile:

```sh
FROM yannvrb56/scrnaseq_r422_seurat43
```

Session info:
```sh
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] fgsea_1.24.0                cluster_2.1.4               gprofiler2_0.2.1           
 [4] edgeR_3.40.2                limma_3.54.2                ComplexHeatmap_2.14.0      
 [7] apeglm_1.20.0               DESeq2_1.38.3               rtracklayer_1.58.0         
[10] ggrepel_0.9.3               XML_3.99-0.13               reprex_2.0.2               
[13] lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0              
[16] purrr_1.0.1                 readr_2.1.4                 tidyr_1.3.0                
[19] tibble_3.2.0                tidyverse_2.0.0             writexl_1.4.2              
[22] readxl_1.4.2                ggpubr_0.6.0                SingleR_2.0.0              
[25] SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
[28] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2           
[31] BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_0.63.0         
[34] openxlsx_4.2.5.2            ggplot2_3.4.1               dplyr_1.1.1                
[37] SeuratObject_4.1.3          Seurat_4.3.0               

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3                scattermore_0.8               coda_0.19-4                  
  [4] bit64_4.0.5                   knitr_1.42                    irlba_2.3.5.1                
  [7] DelayedArray_0.24.0           data.table_1.14.8             KEGGREST_1.38.0              
 [10] RCurl_1.98-1.10               doParallel_1.0.17             generics_0.1.3               
 [13] ScaledMatrix_1.6.0            cowplot_1.1.1                 RSQLite_2.3.0                
 [16] RANN_2.6.1                    future_1.32.0                 bit_4.0.5                    
 [19] tzdb_0.3.0                    spatstat.data_3.0-1           httpuv_1.6.9                 
 [22] xfun_0.37                     celldex_1.8.0                 hms_1.1.2                    
 [25] evaluate_0.20                 promises_1.2.0.1              fansi_1.0.4                  
 [28] restfulr_0.0.15               dbplyr_2.3.1                  igraph_1.4.1                 
 [31] DBI_1.1.3                     geneplotter_1.76.0            htmlwidgets_1.6.1            
 [34] spatstat.geom_3.1-0           ellipsis_0.3.2                backports_1.4.1              
 [37] annotate_1.76.0               deldir_1.0-6                  sparseMatrixStats_1.10.0     
 [40] vctrs_0.6.1                   ROCR_1.0-11                   abind_1.4-5                  
 [43] cachem_1.0.7                  withr_2.5.0                   progressr_0.13.0             
 [46] bdsmatrix_1.3-6               sctransform_0.3.5             GenomicAlignments_1.34.1     
 [49] goftest_1.2-3                 ExperimentHub_2.6.0           lazyeval_0.2.2               
 [52] crayon_1.5.2                  spatstat.explore_3.1-0        pkgconfig_2.0.3              
 [55] nlme_3.1-160                  rlang_1.1.0                   globals_0.16.2               
 [58] lifecycle_1.0.3               miniUI_0.1.1.1                filelock_1.0.2               
 [61] BiocFileCache_2.6.1           rsvd_1.0.5                    AnnotationHub_3.6.0          
 [64] cellranger_1.1.0              polyclip_1.10-4               lmtest_0.9-40                
 [67] Matrix_1.5-1                  carData_3.0-5                 zoo_1.8-11                   
 [70] ggridges_0.5.4                GlobalOptions_0.1.2           png_0.1-8                    
 [73] viridisLite_0.4.1             rjson_0.2.21                  bitops_1.0-7                 
 [76] KernSmooth_2.23-20            Biostrings_2.66.0             blob_1.2.3                   
 [79] DelayedMatrixStats_1.20.0     shape_1.4.6                   parallelly_1.35.0            
 [82] spatstat.random_3.1-4         rstatix_0.7.2                 ggsignif_0.6.4               
 [85] beachmat_2.14.0               scales_1.2.1                  memoise_2.0.1                
 [88] magrittr_2.0.3                plyr_1.8.8                    ica_1.0-3                    
 [91] zlibbioc_1.44.0               compiler_4.2.2                BiocIO_1.8.0                 
 [94] bbmle_1.0.25                  RColorBrewer_1.1-3            clue_0.3-64                  
 [97] fitdistrplus_1.1-8            Rsamtools_2.14.0              cli_3.6.0                    
[100] XVector_0.38.0                listenv_0.9.0                 patchwork_1.1.2              
[103] pbapply_1.7-0                 MASS_7.3-58.1                 tidyselect_1.2.0             
[106] stringi_1.7.12                emdbook_1.3.12                yaml_2.3.7                   
[109] BiocSingular_1.14.0           locfit_1.5-9.7                fastmatch_1.1-3              
[112] tools_4.2.2                   timechange_0.2.0              future.apply_1.10.0          
[115] parallel_4.2.2                circlize_0.4.15               rstudioapi_0.14              
[118] foreach_1.5.2                 gridExtra_2.3                 Rtsne_0.16                   
[121] digest_0.6.31                 BiocManager_1.30.20           shiny_1.7.4                  
[124] Rcpp_1.0.10                   car_3.1-2                     broom_1.0.4                  
[127] BiocVersion_3.16.0            later_1.3.0                   RcppAnnoy_0.0.20             
[130] httr_1.4.5                    AnnotationDbi_1.60.2          colorspace_2.1-0             
[133] fs_1.6.1                      tensor_1.5                    reticulate_1.28              
[136] splines_4.2.2                 uwot_0.1.14                   spatstat.utils_3.0-2         
[139] sp_1.6-0                      plotly_4.10.1                 xtable_1.8-4                 
[142] jsonlite_1.8.4                R6_2.5.1                      pillar_1.8.1                 
[145] htmltools_0.5.4               mime_0.12                     glue_1.6.2                   
[148] fastmap_1.1.1                 BiocParallel_1.32.6           interactiveDisplayBase_1.36.0
[151] codetools_0.2-18              mvtnorm_1.1-3                 utf8_1.2.3                   
[154] lattice_0.20-45               spatstat.sparse_3.0-1         numDeriv_2016.8-1.1          
[157] curl_5.0.0                    leiden_0.4.3                  zip_2.2.2                    
[160] survival_3.4-0                rmarkdown_2.20                munsell_0.5.0                
[163] GetoptLong_1.0.5              GenomeInfoDbData_1.2.9        iterators_1.0.14             
[166] reshape2_1.4.4                gtable_0.3.1
```
