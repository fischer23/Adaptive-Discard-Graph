# Adaptive-Discard-Graph

This repository contains all the code used for the simulations and case study in the paper "An Adaptive-Discard-Graph for online error control" (https://doi.org/10.48550/arXiv.2301.11711). The procedures (ADDIS-Spending_{local}, ADDIS-Graph_{local} and FDR-ADDIS-Graph_{async}) can be found in procedures.R and the code that generates the plots in the Plot_XXX.R files.

In the aforementioned paper, we introduced graphical ADDIS procedures for online multiple testing. The ADDIS-Graph controls the FWER when the p-values are independent. However, it can also be adapted to a local dependence structure and an asynchronous testing setup. In these cases it improves the ADDIS-Spending proposed by Tian & Ramdas (2021). Furthermore, we propose an FDR-ADDIS-Graph, which is superior to the ADDIS* method by Tian and Ramdas (2019). The FDR-ADDIS-Graph controls the FDR under independence of the p-values, however, it can also be extended to an asynchronous setting.
