


The data preparation and precomputations, as well as the analysis of the network, were performed in R, sessionInfo:
```
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
 [6] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] beeswarm_0.4.0        gtable_0.3.5          circlize_0.4.16       shape_1.4.6.1        
 [5] rjson_0.2.21          GlobalOptions_0.1.2   RApiSerialize_0.1.2   tzdb_0.4.0           
 [9] Cairo_1.6-2           vctrs_0.6.5           tools_4.4.0           generics_0.1.3       
[13] stats4_4.4.0          parallel_4.4.0        fansi_1.0.6           cluster_2.1.6        
[17] pkgconfig_2.0.3       RColorBrewer_1.1-3    S4Vectors_0.41.7      RcppParallel_5.1.7   
[21] lifecycle_1.0.4       compiler_4.4.0        farver_2.1.1          textshaping_0.3.7    
[25] munsell_0.5.1         qs_0.26.1             codetools_0.2-20      ComplexHeatmap_2.20.0
[29] clue_0.3-65           vipor_0.4.7           pillar_1.9.0          crayon_1.5.2         
[33] iterators_1.0.14      foreach_1.5.2         tidyselect_1.2.1      ggh4x_0.2.8          
[37] digest_0.6.35         stringi_1.8.3         labeling_0.4.3        grid_4.4.0           
[41] colorspace_2.1-0      cli_3.6.2             magrittr_2.0.3        utf8_1.2.4           
[45] withr_3.0.0           scales_1.3.0          bit64_4.0.5           writexl_1.5.0        
[49] ggbeeswarm_0.7.2      timechange_0.3.0      matrixStats_1.3.0     bit_4.0.5            
[53] ragg_1.3.0            png_0.1-8             GetoptLong_1.0.5      hms_1.1.3            
[57] stringfish_0.16.0     IRanges_2.37.1        doParallel_1.0.17     rlang_1.1.3          
[61] Rcpp_1.0.12           printMat_0.1.0        glue_1.7.0            BiocGenerics_0.50.0  
[65] rstudioapi_0.16.0     vroom_1.6.5           R6_2.5.1              systemfonts_1.0.6   
```



The network computation was performed on a computing cluster. R sessionInfo:
```
R version 4.2.0 (2022-04-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS/LAPACK: /vast/palmer/apps/avx2/software/OpenBLAS/0.3.12-GCC-10.2.0/lib/libopenblas_haswellp-r0.3.12.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] scio_0.9.0      flare_1.7.0.1   igraph_1.3.1    Matrix_1.4-1
 [5] MASS_7.3-57     lattice_0.20-45 QUIC_1.1.1      getopt_1.20.3
 [9] lubridate_1.8.0 forcats_0.5.1   stringr_1.4.0   dplyr_1.1.4
[13] purrr_1.0.2     readr_2.1.5     tidyr_1.2.0     tibble_3.2.1
[17] ggplot2_3.3.5   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.12         pillar_1.9.0        compiler_4.2.0
 [4] tools_4.2.0         qs_0.26.1           lifecycle_1.0.4
 [7] gtable_0.3.0        pkgconfig_2.0.3     rlang_1.1.3
[10] cli_3.6.2           withr_3.0.0         generics_0.1.2
[13] vctrs_0.6.5         hms_1.1.1           grid_4.2.0
[16] tidyselect_1.2.1    glue_1.7.0          RApiSerialize_0.1.2
[19] impute_1.72.3       R6_2.5.1            fansi_1.0.6
[22] projectNPN_0.0.7    glasso_1.11         tzdb_0.4.0
[25] magrittr_2.0.3      scales_1.2.0        ellipsis_0.3.2
[28] colorspace_2.0-3    utf8_1.2.4          stringfish_0.16.0
[31] stringi_1.7.6       RcppParallel_5.1.5  munsell_0.5.0
```

Specifically, used the following implementations for precision matrix estimation:

|  method | function | package | version |
|---------|----------|---------|---------|
| glasso  | glasso   | glasso  | 1.11    |
| QUIC    | QUIC     | QUIC    | 1.1.1   |
| CLIME   | flare    | sugm    | 1.7.0.1 |
| SCIO    | scio     | scio    | 0.9.0   |

The knn imputation was with `{impute}` 1.72.3, installed from BioConductor.



The special-purpose package `{projectNPN}` is available [on Github](https://github.com/cengenproject/projectNPN) to keep all projection and reverse projection functions together with tests.


In addition, the network representations were build using Cytoscape 3.8.2 and R, sessionInfo:

```
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RCy3_2.24.0     wbData_0.9.2    tidygraph_1.3.1 ggraph_2.2.1    igraph_2.0.3   
 [6] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
[11] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1    viridisLite_0.4.2   IRdisplay_1.1       vipor_0.4.7        
 [5] farver_2.1.1        viridis_0.6.5       bitops_1.0-7        fastmap_1.1.1      
 [9] RCurl_1.98-1.14     tweenr_2.0.3        base64url_1.4       XML_3.99-0.16.1    
[13] stringfish_0.16.0   digest_0.6.35       timechange_0.3.0    lifecycle_1.0.4    
[17] magrittr_2.0.3      compiler_4.4.0      rlang_1.1.3         tools_4.4.0        
[21] utf8_1.2.4          labeling_0.4.3      graphlayouts_1.1.1  bit_4.0.5          
[25] RColorBrewer_1.1-3  repr_1.1.7          KernSmooth_2.23-22  pbdZMQ_0.3-11      
[29] withr_3.0.0         BiocGenerics_0.50.0 stats4_4.4.0        grid_4.4.0         
[33] polyclip_1.10-6     fansi_1.0.6         caTools_1.18.2      colorspace_2.1-0   
[37] gtools_3.9.5        scales_1.3.0        MASS_7.3-60.2       cli_3.6.2          
[41] crayon_1.5.2        generics_0.1.3      RcppParallel_5.1.7  rstudioapi_0.16.0  
[45] httr_1.4.7          tzdb_0.4.0          ggbeeswarm_0.7.2    cachem_1.0.8       
[49] ggforce_0.4.2       RApiSerialize_0.1.2 parallel_4.4.0      base64enc_0.1-3    
[53] vctrs_0.6.5         jsonlite_1.8.8      hms_1.1.3           printMat_0.1.0     
[57] bit64_4.0.5         ggrepel_0.9.5       beeswarm_0.4.0      clipr_0.8.0        
[61] glue_1.7.0          RJSONIO_1.3-1.9     stringi_1.8.3       gtable_0.3.5       
[65] munsell_0.5.1       pillar_1.9.0        gplots_3.1.3.1      rappdirs_0.3.3     
[69] htmltools_0.5.8.1   graph_1.82.0        IRkernel_1.3.2      R6_2.5.1           
[73] vroom_1.6.5         evaluate_0.23       backports_1.4.1     qs_0.26.1          
[77] memoise_2.0.1       Rcpp_1.0.12         uuid_1.2-0          gridExtra_2.3      
[81] fs_1.6.4            pkgconfig_2.0.3  
```








