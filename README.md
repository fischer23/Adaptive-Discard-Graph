# Adaptive-Discard-Graph

This repository contains all the code used for the simulations and case study in the paper "An Adaptive-Discard-Graph for online error control" (https://doi.org/10.48550/arXiv.2301.11711). The procedures (ADDIS-Spending_{local}, ADDIS-Graph_{local} and FDR-ADDIS-Graph_{async}) can be found in procedures.R and the code that generates the plots in the Plot_XXX.R files.

In the aforementioned paper, we introduced graphical ADDIS procedures for online multiple testing. The ADDIS-Graph and FDR-ADDIS-Graph control the FWER and FDR, respectively, when the null p-values are independent from each other and the non-nulls. However, both procedures can be extended to a local dependence structure and an asynchronous testing setup, while the FDR-ADDIS-Graph provides only marginal FDR control under local dependence. 
