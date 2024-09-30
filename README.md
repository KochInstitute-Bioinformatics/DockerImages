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

Run docker build command to build your image:

```sh
docker build . -t yannvrb56/bulkrnaseq
```
# DockerImages

## Bulk RNA-seq

To pull the image from the command line:

```sh
docker pull yannvrb56/bulkrnaseq
```

To use as a base image in a new Dockerfile:

```sh
FROM yannvrb56/bulkrnaseq
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

## Single cell RNA-seq and single cell ATAC-seq (Seurat5.0 anf Signac)

yannvrb56/r441seurat5signac

## Single cell RNA-seq (Seurat4.3)

yannvrb56/scrnaseq_r422_seurat43
