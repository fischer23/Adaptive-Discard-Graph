# Generates results for Table 1 of the paper
# "ADDIS-Graphs for online error control with application to platform trials"

rm(list = ls())
library(ggplot2)
library(MASS)
library(patchwork)
library(Matrix)
library(mvtnorm)

### load the procedures
source("Graph_Procedures.R")

names <- c(
  "Dexamethasone", "Lopinavir-ritonavir", "Hydroxychloroquine", "Azithromycin",
  "Tocilizumab", "Convalescent plasma", "Casirivimab-Imdevimab", "Aspirin",
  "Colchicine", "Baricitinib", "High-dose steroids", "Empagliflozin"
)

lags <- c(0, 1, 2, 3, 4, 5, 3, 3, 3, 3, 1, 2)

p <- c(0.0003, 0.58, 0.1, 0.99, 0.007, 0.34, 0.001, 0.35, 0.63, 0.026, 0.0012, 0.64)

n <- 12

### We run the experiment for different gamma
n_q <- 3

alpha_ind_unc <- matrix(0, ncol = n_q, nrow = n)
rej_unc <- matrix(0, ncol = n_q, nrow = n)
n_rej_unc <- rep(0, n_q)
future_level_unc <- rep(0, n_q)
alpha_ind_Spending <- matrix(0, ncol = n_q, nrow = n)
rej_Spending <- matrix(0, ncol = n_q, nrow = n)
n_rej_Spending <- rep(0, n_q)
future_level_Spending <- rep(0, n_q)
alpha_ind_Graph <- matrix(0, ncol = n_q, nrow = n)
rej_Graph <- matrix(0, ncol = n_q, nrow = n)
n_rej_Graph <- rep(0, n_q)
future_level_Graph <- rep(0, n_q)

count <- 1

for (q in c(0.6, 0.7, 0.8)) {
  ### Initialise Hyperparameters
  # gamma=0.5^(1:n)      #2.062 is the approximated value such that the series equals 1
  gamma <- q^(1:n) * (1 - q) / q
  tau <- 0.8
  lambda <- 0.3
  lambda_corr <- 0.5
  w <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  w[w == 0] <- 1
  w <- matrix(gamma[w], n, n)
  w[upper.tri(w) == 0] <- 0
  alpha <- 0.05

  ### Calculation d_j and weights for Adaptive-Graph_{corr}

  # d_j for ADDIS-Graph_{local-u}
  k_lag <- seq(1, n, 1) - lags
  d_j <- rep((n + 1), n)
  for (k in 1:n) {
    for (i in k:n) {
      if (k_lag[i] > k) {
        d_j[k] <- i
        break
      }
    }
  }

  ## Uncorrrected Testing
  alpha_ind_unc[, count] <- rep(alpha, n)
  rej_unc[, count] <- alpha_ind_unc[, count] >= p
  n_rej_unc[count] <- sum(rej_unc[, count])
  future_level_unc[count] <- Inf

  ## ADDIS-Spending
  alpha_ind_Spending[, count] <- ADDIS_Spending(alpha, gamma, tau, lambda, lags, p, n)
  rej_Spending[, count] <- alpha_ind_Spending[, count] >= p
  n_rej_Spending[count] <- sum(rej_Spending[, count])
  future_level_Spending[count] <- alpha * sum(q^((sum((p > lambda & p <= tau)) + 1):10000) * (1 - q) / q)

  ## ADDIS-Graph_{local-u}
  alpha_ind_Graph[, count] <- ADDIS_Graph_imp(alpha, gamma, tau, lambda, lags, d_j, p, n)
  rej_Graph[, count] <- alpha_ind_Graph[, count] >= p
  n_rej_Graph[count] <- sum(rej_Graph[, count])
  future_level_Graph[count] <- alpha - (sum(alpha_ind_Graph[which((p <= lambda | p > tau) == 0), count])) / (tau - lambda)

  count <- count + 1
}
