### This R-file generates the plots of the paper
# "ADDIS-Graphs for online error control with application to platform trials"

library(ggplot2)
library(patchwork)

# Load plot functions

source("sim_functions.R")

### Figure 8

p1 <- plot_generator_FWER_batchsize("results/FWER_gamma_logq.rds", "Spending")

p2 <- plot_generator_FWER_batchsize("results/FWER_gamma_q_1_6.rds", "Spending")

p3 <- plot_generator_FWER_batchsize("results/FWER_gamma_q_2.rds", "Spending")

combined_plot <- p1 / p2 / p3

ggsave("results/Figure8.pdf", combined_plot, width = 7.64, height = 10.5)

### Figure S.1

p1 <- plot_generator_FWER_batchsize("results/FWER_gamma_logq.rds", "Spending_c")

p2 <- plot_generator_FWER_batchsize("results/FWER_gamma_q_1_6.rds", "Spending_c")

p3 <- plot_generator_FWER_batchsize("results/FWER_gamma_q_2.rds", "Spending_c")

combined_plot <- p1 / p2 / p3

ggsave("results/FigureS1.pdf", combined_plot, width = 7.64, height = 10.5)

### Figure S.2

p1 <- plot_generator_FWER_batchsize("results/FWER_mu_N_0.rds", "Spending")

p2 <- plot_generator_FWER_batchsize("results/FWER_mu_N_1.rds", "Spending")

p3 <- plot_generator_FWER_batchsize("results/FWER_mu_N_2.rds", "Spending")

combined_plot <- p1 / p2 / p3

ggsave("results/FigureS2.pdf", combined_plot, width = 7.64, height = 10.5)

### Figure S.3

p1 <- plot_generator_FWER_batchsize("results/FWER_rho_2.rds", "Spending")

p2 <- plot_generator_FWER_batchsize("results/FWER_rho_8.rds", "Spending")

combined_plot <- p1 / p2

ggsave("results/FigureS3.pdf", combined_plot, width = 7.64, height = 7)

### Figure S.4

p1 <- plot_generator_FWER_smallest_batch("results/FWER_incr_batches.rds")

p2 <- plot_generator_FWER_smallest_batch("results/FWER_decr_batches.rds")

combined_plot <- p1 / p2

ggsave("results/FigureS4.pdf", combined_plot, width = 7.64, height = 7)

### Figure S.5

p1 <- plot_generator_FWER_corr("results/FWER_corr_batch.rds", "batch")

p2 <- plot_generator_FWER_corr("results/FWER_corr_rho.rds", "rho")

p3 <- plot_generator_FWER_corr("results/FWER_corr_mu_N.rds", "mu_N")

combined_plot <- p1 / p2 / p3

ggsave("results/FigureS5.pdf", combined_plot, width = 7.64, height = 10.5)

### Figure S.6

p1 <- plot_generator_FDR_async("results/FDR_gamma_logq.rds")

p2 <- plot_generator_FDR_async("results/FDR_gamma_q_1_6.rds")

p3 <- plot_generator_FDR_async("results/FDR_gamma_q_2.rds")

combined_plot <- p1 / p2 / p3

ggsave("results/FigureS6.pdf", combined_plot, width = 7.64, height = 10.5)
