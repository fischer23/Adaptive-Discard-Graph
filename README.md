# Adaptive-Discard-Graph

This repository contains all the code used for the simulations and case study in the paper "ADDIS-Graphs for online
error control with application to platform trials" (https://arxiv.org/abs/2301.11711). In the aforementioned paper we introduce a graphical Adaptive-Discard procedure for online error control. Our ADDIS-Graph is particularly well suited for online settings with conflicts such as encountered in platform trials. The ADDIS-Graph uniformly improves the ADDIS-Spending by Tian and Ramdas (2021) when conflicts are present.

Files:

Graph_Procedures.R:  Contains functions that implement all procedures of our paper.

sim_functions.R:     Contains functions that implement the simulation settings considered in our paper.

generate_sim_results.R: Uses the functions in sim_functions.R to generate all simulations results and saves them in the results folder.

generate_plots.R: Uses the simulation results in the results folder to create all plots of our paper.

real_data_RECOVERY.R: Applies our ADDIS-Graph to the real data of the RECOVERY trial. The used p-values were extracted from the RECOVERY trial website https://www.recoverytrial.net/results and the local dependence structure by a publication of the data monitoring commitee of the RECOVERY trial (Sandercock et al. 2022).

results folder: Contains all results of the paper.

R version 4.5.0 (2025-04-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8
[4] LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] mvtnorm_1.3-3   Matrix_1.7-3    reshape_0.8.9   MASS_7.3-65     patchwork_1.3.0 ggplot2_3.5.2  

loaded via a namespace (and not attached):
 [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        generics_0.1.4     glue_1.8.0        
 [6] plyr_1.8.9         scales_1.4.0       grid_4.5.0         tibble_3.2.1       lifecycle_1.0.4   
[11] compiler_4.5.0     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.0.14        pkgconfig_2.0.3   
[16] rstudioapi_0.17.1  farver_2.1.2       lattice_0.22-7     R6_2.6.1           tidyselect_1.2.1  
[21] pillar_1.10.2      magrittr_2.0.3     tools_4.5.0        withr_3.0.2        gtable_0.3.6  




