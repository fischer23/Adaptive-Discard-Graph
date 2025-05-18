# Generates plots for Figure 8 of the paper
# "ADDIS-Graphs for onlineerror control with application to platform trials"

rm(list = ls())
library(ggplot2)
library(MASS)
library(patchwork)

### load the procedures
source("Graph_Procedures.R")

### Gaussian testing problem for exact p-values
m <- 1000 # Number of Trials
n <- 100 # Number of Hypotheses per Trial
mu_A <- 3 # Strength of the alternative
mu_N <- -0.5 # Conservativeness of null p-values (<0 for conservative null p-values)
# pi_A is defined in the loop below

### Initialise Hyperparameters
seq_1_n <- seq(1, n, 1)
gamma <- (1 / ((seq_1_n + 1) * log(seq_1_n + 1)^2)) / 2.06227 # 2.06227 is the approximated value such that the series equals 1
alpha <- 0.2
tau <- 0.8
lambda <- 0.16

# Correlation within batch
corr <- 0.5

### Set seed to make the results reproducible
set.seed(12345)

# Set the different batch sizes
batch_sizes <- c(1, 5, 10, 20)

### Predefine vectors for FWER and power of the different procedures
FWER_Spending <- matrix(0, 9, length(batch_sizes))
power_Spending <- matrix(0, 9, length(batch_sizes))
FWER_Graph <- matrix(0, 9, length(batch_sizes))
power_Graph <- matrix(0, 9, length(batch_sizes))
b <- 1 # Counter
for (batch_size in batch_sizes) {
  ### Parameters for a local dependence structure given by batches
  batch_number <- ceiling(n / batch_size)
  lags <- rep(seq(0, (batch_size - 1), 1), batch_number)
  lags <- lags[1:n]
  sigma <- matrix(corr, batch_size, batch_size) + diag((1 - corr), batch_size)
  mu <- rep(0, batch_size)

  # Calculate d_j
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


  ### Generate p-values and compute FWER and power for the desired procedures
  for (l in 1:9) {
    ## Fast way for p-values in batches.
    pi_A <- l / 10
    p <- matrix(, nrow = n, ncol = m)
    hypo <- matrix(, nrow = n, ncol = m)
    for (j in 1:m) {
      hypo[, j] <- rbinom(n, 1, pi_A)
      X <- rep(0, batch_number * batch_size)
      for (k in 1:batch_number) {
        X[((k - 1) * batch_size + 1):(k * batch_size)] <- mvrnorm(1, mu, sigma)
      }
      X <- X[1:n]
      Z <- mu_N * (hypo[, j] - 1) * (-1) + mu_A * hypo[, j] + X
      p[, j] <- pnorm(-Z)
    }

    ## ADDIS-Spending
    V <- rep(0, m) # Indicates, whether there was atleast one type 1 error in a trial
    power <- rep(0, m) # Power within each trial

    for (j in 1:m) {
      alpha_ind <- closed_ADDIS_Spending(alpha, gamma, tau, lambda, lags, p[, j], n)
      hypo_est <- alpha_ind >= p[, j]
      V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
      D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
      power[j] <- sum(D) / sum(hypo[, j])
    }
    FWER_Spending[l, b] <- mean(V)
    power_Spending[l, b] <- mean(power, na.rm = TRUE)


    ## ADDIS-Graph_{local-u}
    V <- rep(0, m) # Indicates, whether there was atleast one type 1 error in a trial
    power <- rep(0, m) # Power within each trial

    for (j in 1:m) {
      alpha_ind <- ADDIS_Graph_imp(alpha, gamma, tau, lambda, lags, d_j, p[, j], n)
      hypo_est <- alpha_ind >= p[, j]
      V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
      D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
      power[j] <- sum(D) / sum(hypo[, j])
    }
    FWER_Graph[l, b] <- mean(V)
    power_Graph[l, b] <- mean(power, na.rm = TRUE)
  }
  b <- b + 1
}
### Create Plot for ADDIS-Spending

results_df <- data.frame(seq(0.1, 0.9, 0.1), power_Spending, FWER_Spending)

p_spending <- ggplot(results_df, aes(seq.0.1..0.9..0.1.)) +
  geom_line(aes(y = X1, linetype = "1", colour = "2")) +
  geom_point(aes(y = X1, colour = "2")) +
  geom_line(aes(y = X1.1, linetype = "1", colour = "2")) +
  geom_point(aes(y = X1.1, colour = "2")) +
  geom_line(aes(y = X2, linetype = "2", colour = "2")) +
  geom_point(aes(y = X2, colour = "2")) +
  geom_line(aes(y = X2.1, linetype = "2", colour = "2")) +
  geom_point(aes(y = X2.1, colour = "2")) +
  geom_line(aes(y = X3, linetype = "3", colour = "2")) +
  geom_point(aes(y = X3, colour = "2")) +
  geom_line(aes(y = X3.1, linetype = "3", colour = "2")) +
  geom_point(aes(y = X3.1, colour = "2")) +
  geom_line(aes(y = X4, linetype = "4", colour = "2")) +
  geom_point(aes(y = X4, colour = "2")) +
  geom_line(aes(y = X4.1, linetype = "4", colour = "2")) +
  geom_point(aes(y = X4.1, colour = "2")) +
  geom_hline(yintercept = alpha) +
  scale_colour_manual(guide = "none", values = c("2" = "#f84f4f")) +
  scale_linetype_manual(
    name = "Batch-size", values = c("1" = "solid", "2" = "longdash", "3" = "dashed", "4" = "dotted"),
    labels = c("1", "5", "10", "20")
  ) +
  xlab(expression(pi[A])) +
  ylab("FWER / Power") +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.2), limits = c(0.1, 0.9), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  theme(
    panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  ggtitle(expression("Closed ADDIS-Spending"[local]))


### Create Plot for ADDIS-Graph
results_df_graph <- data.frame(seq(0.1, 0.9, 0.1), power_Graph, FWER_Graph)

p_graph <- ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) +
  geom_line(aes(y = X1, linetype = "1", colour = "4")) +
  geom_point(aes(y = X1, colour = "4", shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "1", colour = "4")) +
  geom_point(aes(y = X1.1), colour = "4", shape = "4") +
  geom_line(aes(y = X2, linetype = "2", colour = "4")) +
  geom_point(aes(y = X2, colour = "4", shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2", colour = "4")) +
  geom_point(aes(y = X2.1, colour = "4", shape = "4")) +
  geom_line(aes(y = X3, linetype = "3", colour = "4")) +
  geom_point(aes(y = X3, colour = "4", shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "3", colour = "4")) +
  geom_point(aes(y = X3.1, colour = "4", shape = "4")) +
  geom_line(aes(y = X4, linetype = "4", colour = "4")) +
  geom_point(aes(y = X4, colour = "4", shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "4", colour = "4")) +
  geom_point(aes(y = X4.1, colour = "4", shape = "4")) +
  geom_hline(yintercept = alpha) +
  scale_shape_manual(guide = "none", values = c("4" = 17)) +
  scale_colour_manual(guide = "none", values = c("4" = "cornflowerblue")) +
  scale_linetype_manual(
    name = "Batch-size", values = c("1" = "solid", "2" = "longdash", "3" = "dashed", "4" = "dotted"),
    labels = c("1", "5", "10", "20")
  ) +
  xlab(expression(pi[A])) +
  ylab("FWER / Power") +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.2), limits = c(0.1, 0.9), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  theme(
    panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  ggtitle(expression("ADDIS-Graph"[conf - u]))

combined <- p_spending + p_graph & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FWER_closed_ADDIS_Graph_local-u_logq.pdf", width = 7.64, height = 3.59)


############################# Results for gamma_i=6/(pi^2*i^2)
gamma <- 6 / (pi^2 * (seq_1_n^2))

### Predefine vectors for FWER and power of the different procedures
FWER_Spending <- matrix(0, 9, length(batch_sizes))
power_Spending <- matrix(0, 9, length(batch_sizes))
FWER_Graph <- matrix(0, 9, length(batch_sizes))
power_Graph <- matrix(0, 9, length(batch_sizes))
b <- 1 # Counter
for (batch_size in batch_sizes) {
  ### Parameters for a local dependence structure given by batches
  batch_number <- ceiling(n / batch_size)
  lags <- rep(seq(0, (batch_size - 1), 1), batch_number)
  lags <- lags[1:n]
  sigma <- matrix(corr, batch_size, batch_size) + diag((1 - corr), batch_size)
  mu <- rep(0, batch_size)

  # Calculate d_j
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


  ### Generate p-values and compute FWER and power for the desired procedures
  for (l in 1:9) {
    ## Fast way for p-values in batches.
    pi_A <- l / 10
    p <- matrix(, nrow = n, ncol = m)
    hypo <- matrix(, nrow = n, ncol = m)
    for (j in 1:m) {
      hypo[, j] <- rbinom(n, 1, pi_A)
      X <- rep(0, batch_number * batch_size)
      for (k in 1:batch_number) {
        X[((k - 1) * batch_size + 1):(k * batch_size)] <- mvrnorm(1, mu, sigma)
      }
      X <- X[1:n]
      Z <- mu_N * (hypo[, j] - 1) * (-1) + mu_A * hypo[, j] + X
      p[, j] <- pnorm(-Z)
    }

    ## ADDIS-Spending
    V <- rep(0, m) # Indicates, whether there was atleast one type 1 error in a trial
    power <- rep(0, m) # Power within each trial

    for (j in 1:m) {
      alpha_ind <- closed_ADDIS_Spending(alpha, gamma, tau, lambda, lags, p[, j], n)
      hypo_est <- alpha_ind >= p[, j]
      V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
      D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
      power[j] <- sum(D) / sum(hypo[, j])
    }
    FWER_Spending[l, b] <- mean(V)
    power_Spending[l, b] <- mean(power, na.rm = TRUE)


    ## ADDIS-Graph_{local-u}
    V <- rep(0, m) # Indicates, whether there was atleast one type 1 error in a trial
    power <- rep(0, m) # Power within each trial

    for (j in 1:m) {
      alpha_ind <- ADDIS_Graph_imp(alpha, gamma, tau, lambda, lags, d_j, p[, j], n)
      hypo_est <- alpha_ind >= p[, j]
      V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
      D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
      power[j] <- sum(D) / sum(hypo[, j])
    }
    FWER_Graph[l, b] <- mean(V)
    power_Graph[l, b] <- mean(power, na.rm = TRUE)
  }
  b <- b + 1
}
### Create Plot for ADDIS-Spending

results_df <- data.frame(seq(0.1, 0.9, 0.1), power_Spending, FWER_Spending)

p_spending <- ggplot(results_df, aes(seq.0.1..0.9..0.1.)) +
  geom_line(aes(y = X1, linetype = "1", colour = "2")) +
  geom_point(aes(y = X1, colour = "2")) +
  geom_line(aes(y = X1.1, linetype = "1", colour = "2")) +
  geom_point(aes(y = X1.1, colour = "2")) +
  geom_line(aes(y = X2, linetype = "2", colour = "2")) +
  geom_point(aes(y = X2, colour = "2")) +
  geom_line(aes(y = X2.1, linetype = "2", colour = "2")) +
  geom_point(aes(y = X2.1, colour = "2")) +
  geom_line(aes(y = X3, linetype = "3", colour = "2")) +
  geom_point(aes(y = X3, colour = "2")) +
  geom_line(aes(y = X3.1, linetype = "3", colour = "2")) +
  geom_point(aes(y = X3.1, colour = "2")) +
  geom_line(aes(y = X4, linetype = "4", colour = "2")) +
  geom_point(aes(y = X4, colour = "2")) +
  geom_line(aes(y = X4.1, linetype = "4", colour = "2")) +
  geom_point(aes(y = X4.1, colour = "2")) +
  geom_hline(yintercept = alpha) +
  scale_colour_manual(guide = "none", values = c("2" = "#f84f4f")) +
  scale_linetype_manual(
    name = "Batch-size", values = c("1" = "solid", "2" = "longdash", "3" = "dashed", "4" = "dotted"),
    labels = c("1", "5", "10", "20")
  ) +
  xlab(expression(pi[A])) +
  ylab("FWER / Power") +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.2), limits = c(0.1, 0.9), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  theme(
    panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  ggtitle(expression("Closed ADDIS-Spending"[local]))


### Create Plot for ADDIS-Graph
results_df_graph <- data.frame(seq(0.1, 0.9, 0.1), power_Graph, FWER_Graph)

p_graph <- ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) +
  geom_line(aes(y = X1, linetype = "1", colour = "4")) +
  geom_point(aes(y = X1, colour = "4", shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "1", colour = "4")) +
  geom_point(aes(y = X1.1), colour = "4", shape = "4") +
  geom_line(aes(y = X2, linetype = "2", colour = "4")) +
  geom_point(aes(y = X2, colour = "4", shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2", colour = "4")) +
  geom_point(aes(y = X2.1, colour = "4", shape = "4")) +
  geom_line(aes(y = X3, linetype = "3", colour = "4")) +
  geom_point(aes(y = X3, colour = "4", shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "3", colour = "4")) +
  geom_point(aes(y = X3.1, colour = "4", shape = "4")) +
  geom_line(aes(y = X4, linetype = "4", colour = "4")) +
  geom_point(aes(y = X4, colour = "4", shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "4", colour = "4")) +
  geom_point(aes(y = X4.1, colour = "4", shape = "4")) +
  geom_hline(yintercept = alpha) +
  scale_shape_manual(guide = "none", values = c("4" = 17)) +
  scale_colour_manual(guide = "none", values = c("4" = "cornflowerblue")) +
  scale_linetype_manual(
    name = "Batch-size", values = c("1" = "solid", "2" = "longdash", "3" = "dashed", "4" = "dotted"),
    labels = c("1", "5", "10", "20")
  ) +
  xlab(expression(pi[A])) +
  ylab("FWER / Power") +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.2), limits = c(0.1, 0.9), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  theme(
    panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  ggtitle(expression("ADDIS-Graph"[conf - u]))

combined <- p_spending + p_graph & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FWER_closed_ADDIS_Graph_local-u_q^2.pdf", width = 7.64, height = 3.59)


############################# Results for gamma_i\propto 1/(pi^1.6)
gamma <- 1 / (2.28577 * seq_1_n^1.6) # 2.28577 approximated value such that sum of gamma equals one

### Predefine vectors for FWER and power of the different procedures
FWER_Spending <- matrix(0, 9, length(batch_sizes))
power_Spending <- matrix(0, 9, length(batch_sizes))
FWER_Graph <- matrix(0, 9, length(batch_sizes))
power_Graph <- matrix(0, 9, length(batch_sizes))
b <- 1 # Counter
for (batch_size in batch_sizes) {
  ### Parameters for a local dependence structure given by batches
  batch_number <- ceiling(n / batch_size)
  lags <- rep(seq(0, (batch_size - 1), 1), batch_number)
  lags <- lags[1:n]
  sigma <- matrix(corr, batch_size, batch_size) + diag((1 - corr), batch_size)
  mu <- rep(0, batch_size)

  # Calculate d_j
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


  ### Generate p-values and compute FWER and power for the desired procedures
  for (l in 1:9) {
    ## Fast way for p-values in batches.
    pi_A <- l / 10
    p <- matrix(, nrow = n, ncol = m)
    hypo <- matrix(, nrow = n, ncol = m)
    for (j in 1:m) {
      hypo[, j] <- rbinom(n, 1, pi_A)
      X <- rep(0, batch_number * batch_size)
      for (k in 1:batch_number) {
        X[((k - 1) * batch_size + 1):(k * batch_size)] <- mvrnorm(1, mu, sigma)
      }
      X <- X[1:n]
      Z <- mu_N * (hypo[, j] - 1) * (-1) + mu_A * hypo[, j] + X
      p[, j] <- pnorm(-Z)
    }

    ## ADDIS-Spending
    V <- rep(0, m) # Indicates, whether there was atleast one type 1 error in a trial
    power <- rep(0, m) # Power within each trial

    for (j in 1:m) {
      alpha_ind <- closed_ADDIS_Spending(alpha, gamma, tau, lambda, lags, p[, j], n)
      hypo_est <- alpha_ind >= p[, j]
      V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
      D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
      power[j] <- sum(D) / sum(hypo[, j])
    }
    FWER_Spending[l, b] <- mean(V)
    power_Spending[l, b] <- mean(power, na.rm = TRUE)


    ## ADDIS-Graph_{local-u}
    V <- rep(0, m) # Indicates, whether there was atleast one type 1 error in a trial
    power <- rep(0, m) # Power within each trial

    for (j in 1:m) {
      alpha_ind <- ADDIS_Graph_imp(alpha, gamma, tau, lambda, lags, d_j, p[, j], n)
      hypo_est <- alpha_ind >= p[, j]
      V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
      D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
      power[j] <- sum(D) / sum(hypo[, j])
    }
    FWER_Graph[l, b] <- mean(V)
    power_Graph[l, b] <- mean(power, na.rm = TRUE)
  }
  b <- b + 1
}
### Create Plot for ADDIS-Spending

results_df <- data.frame(seq(0.1, 0.9, 0.1), power_Spending, FWER_Spending)

p_spending <- ggplot(results_df, aes(seq.0.1..0.9..0.1.)) +
  geom_line(aes(y = X1, linetype = "1", colour = "2")) +
  geom_point(aes(y = X1, colour = "2")) +
  geom_line(aes(y = X1.1, linetype = "1", colour = "2")) +
  geom_point(aes(y = X1.1, colour = "2")) +
  geom_line(aes(y = X2, linetype = "2", colour = "2")) +
  geom_point(aes(y = X2, colour = "2")) +
  geom_line(aes(y = X2.1, linetype = "2", colour = "2")) +
  geom_point(aes(y = X2.1, colour = "2")) +
  geom_line(aes(y = X3, linetype = "3", colour = "2")) +
  geom_point(aes(y = X3, colour = "2")) +
  geom_line(aes(y = X3.1, linetype = "3", colour = "2")) +
  geom_point(aes(y = X3.1, colour = "2")) +
  geom_line(aes(y = X4, linetype = "4", colour = "2")) +
  geom_point(aes(y = X4, colour = "2")) +
  geom_line(aes(y = X4.1, linetype = "4", colour = "2")) +
  geom_point(aes(y = X4.1, colour = "2")) +
  geom_hline(yintercept = alpha) +
  scale_colour_manual(guide = "none", values = c("2" = "#f84f4f")) +
  scale_linetype_manual(
    name = "Batch-size", values = c("1" = "solid", "2" = "longdash", "3" = "dashed", "4" = "dotted"),
    labels = c("1", "5", "10", "20")
  ) +
  xlab(expression(pi[A])) +
  ylab("FWER / Power") +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.2), limits = c(0.1, 0.9), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  theme(
    panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  ggtitle(expression("Closed ADDIS-Spending"[local]))


### Create Plot for ADDIS-Graph
results_df_graph <- data.frame(seq(0.1, 0.9, 0.1), power_Graph, FWER_Graph)

p_graph <- ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) +
  geom_line(aes(y = X1, linetype = "1", colour = "4")) +
  geom_point(aes(y = X1, colour = "4", shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "1", colour = "4")) +
  geom_point(aes(y = X1.1), colour = "4", shape = "4") +
  geom_line(aes(y = X2, linetype = "2", colour = "4")) +
  geom_point(aes(y = X2, colour = "4", shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2", colour = "4")) +
  geom_point(aes(y = X2.1, colour = "4", shape = "4")) +
  geom_line(aes(y = X3, linetype = "3", colour = "4")) +
  geom_point(aes(y = X3, colour = "4", shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "3", colour = "4")) +
  geom_point(aes(y = X3.1, colour = "4", shape = "4")) +
  geom_line(aes(y = X4, linetype = "4", colour = "4")) +
  geom_point(aes(y = X4, colour = "4", shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "4", colour = "4")) +
  geom_point(aes(y = X4.1, colour = "4", shape = "4")) +
  geom_hline(yintercept = alpha) +
  scale_shape_manual(guide = "none", values = c("4" = 17)) +
  scale_colour_manual(guide = "none", values = c("4" = "cornflowerblue")) +
  scale_linetype_manual(
    name = "Batch-size", values = c("1" = "solid", "2" = "longdash", "3" = "dashed", "4" = "dotted"),
    labels = c("1", "5", "10", "20")
  ) +
  xlab(expression(pi[A])) +
  ylab("FWER / Power") +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.2), limits = c(0.1, 0.9), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  theme(
    panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  ggtitle(expression("ADDIS-Graph"[conf - u]))

combined <- p_spending + p_graph & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FWER_closed_ADDIS_Graph_local-u_q^1,6.pdf", width = 7.64, height = 3.59)
