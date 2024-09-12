# Adaptive-Discard-Graph

This repository contains all the code used for the simulations and case study in the paper "ADDIS-Graphs for online
error control with application to platform trials" (https://arxiv.org/abs/2301.11711). The proposed online procedures can be found in Graph_procedures.R and the code that generates the plots in the Plot_XXX.R files. The real data experiment can be performed with real_data_RECOVERY. The data used in the real data experiment was derived from the RECOVERY trial website https://www.recoverytrial.net/results. 

In the aforementioned paper, we introduced graphical ADDIS procedures for online multiple testing that are particularly useful for platform trials. The ADDIS-Graph and FDR-ADDIS-Graph control the FWER and FDR, respectively, when the null p-values are independent from each other and the non-nulls. However, both procedures can be extended to a local dependence structure and an asynchronous testing setup, whereby the FDR-ADDIS-Graph provides only marginal FDR control under local dependence.
