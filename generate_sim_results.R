library(MASS)
library(patchwork)
library(reshape)
library(Matrix)
library(mvtnorm)

source("sim_functions.R")
source("Graph_Procedures.R")

# Generate sim results used in Figure 8 and S.1

data_generator_FWER_batchsize(
  gamma = (1 / (((1:100) + 1) * log((1:100) + 1)^2)) / 2.06227, mu_N = -0.5,
  corr = 0.5, filename = "results/FWER_gamma_logq.rds"
)

data_generator_FWER_batchsize(
  gamma = (1 / (2.28577 * (1:100)^1.6)), mu_N = -0.5,
  corr = 0.5, filename = "results/FWER_gamma_q_1_6.rds"
)

data_generator_FWER_batchsize(
  gamma = (6 / (pi^2 * (1:100)^2)), mu_N = -0.5,
  corr = 0.5, filename = "results/FWER_gamma_q_2.rds"
)

# Generate sim results used in Figure S.2

data_generator_FWER_batchsize(
  gamma = (6 / (pi^2 * (1:100)^2)), mu_N = 0,
  corr = 0.5, filename = "results/FWER_mu_N_0.rds"
)

data_generator_FWER_batchsize(
  gamma = (6 / (pi^2 * (1:100)^2)), mu_N = -1,
  corr = 0.5, filename = "results/FWER_mu_N_1.rds"
)

data_generator_FWER_batchsize(
  gamma = (6 / (pi^2 * (1:100)^2)), mu_N = -2,
  corr = 0.5, filename = "results/FWER_mu_N_2.rds"
)

# Generate sim results used in Figure S.3

data_generator_FWER_batchsize(
  gamma = (6 / (pi^2 * (1:100)^2)), mu_N = -0.5,
  corr = 0.2, filename = "results/FWER_rho_2.rds"
)

data_generator_FWER_batchsize(
  gamma = (6 / (pi^2 * (1:100)^2)), mu_N = -0.5,
  corr = 0.8, filename = "results/FWER_rho_8.rds"
)

# Generate sim results used in Figure S.4

data_generator_FWER_smallest_batch(smallest_batches = c(1, 2, 5, 10), direction = "incr",
                                   filename = "results/FWER_incr_batches.rds")

data_generator_FWER_smallest_batch(smallest_batches = c(1, 2, 5, 10), direction = "decr",
                                   filename = "results/FWER_decr_batches.rds")

# Generate sim results used in Figure S.5

data_generator_FWER_corr(var_ind = "mu_N", mu_Ns = c(0, -0.5, -1, -2), 
                         filename = "results/FWER_corr_mu_N.rds")

data_generator_FWER_corr(var_ind = "rho", rhos = c(0.3, 0.5, 0.7, 0.9), 
                         filename = "results/FWER_corr_rho.rds")

data_generator_FWER_corr(var_ind = "batch", batch_sizes = c(1, 5, 10, 20), 
                         filename = "results/FWER_corr_batch.rds")

# Generate sim results used in Figure S.6

data_generator_FDR_async(
  gamma = (1 / (((1:100) + 1) * log((1:100) + 1)^2)) / 2.06227, 
  filename = "results/FDR_gamma_logq.rds"
)

data_generator_FDR_async(
  gamma = (1 / (2.28577 * (1:100)^1.6)), 
  filename = "results/FDR_gamma_q_1_6.rds"
)

data_generator_FDR_async(
  gamma = (6 / (pi^2 * (1:100)^2)), 
  filename = "results/FDR_gamma_q_2.rds"
)
