# These functions can generate the simulation results for the different settings in the paper
# "ADDIS-Graphs for online error control with application to platform trials"


### Function to generate the data in the FWER case for different batch sizes

# Input:
# gamma:    initial weights (nonnegative n-dim. vector)
# mu_N:     conservativeness of null p-values (real number <=0)
# corr:     correlation within batch (real number in [0,1))
# filename: Output filename

# Output: Saves simulation results (power and FWER) of the ADDIS-Spending, closed ADDIS-Spending and
#         ADDIS-Graph for the different batch_sizes (1, 5, 10, 20) in file "filename"

data_generator_FWER_batchsize <- function(gamma, mu_N, corr, filename) {
  ### Simulation parameters
  m <- 1000 # Number of Trials
  n <- 100 # Number of hypotheses
  mu_A <- 3 # Strength of the alternative
  # pi_A is defined in the loop below

  ### Initialise Hyperparameters
  alpha <- 0.2
  tau <- 0.8
  lambda <- 0.16

  ### Set seed to make the results reproducible
  set.seed(12345)

  # Set the different batch sizes
  batch_sizes <- c(1, 5, 10, 20)

  ### Predefine vectors for FWER and power of the different procedures
  FWER_Spending <- matrix(0, 9, length(batch_sizes))
  power_Spending <- matrix(0, 9, length(batch_sizes))
  FWER_Spending_c <- matrix(0, 9, length(batch_sizes))
  power_Spending_c <- matrix(0, 9, length(batch_sizes))
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

    # Calculate d_j for ADDIS-Graph
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
      # Fast way for p-values in batches.
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
        Z <- mu_N * (1 - hypo[, j]) + mu_A * hypo[, j] + X
        p[, j] <- pnorm(-Z)
      }

      ## ADDIS-Spending
      V <- rep(0, m) # Indicates, whether there was at least one type 1 error in a trial
      power <- rep(0, m) # Power within each trial

      for (j in 1:m) {
        alpha_ind <- ADDIS_Spending(alpha, gamma, tau, lambda, lags, p[, j], n)
        hypo_est <- alpha_ind >= p[, j]
        V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
        D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
        power[j] <- sum(D) / sum(hypo[, j])
      }
      FWER_Spending[l, b] <- mean(V)
      power_Spending[l, b] <- mean(power, na.rm = TRUE)

      ## Closed ADDIS-Spending
      V <- rep(0, m) # Indicates, whether there was at least one type 1 error in a trial
      power <- rep(0, m) # Power within each trial

      for (j in 1:m) {
        alpha_ind <- closed_ADDIS_Spending(alpha, gamma, tau, lambda, lags, p[, j], n)
        hypo_est <- alpha_ind >= p[, j]
        V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
        D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
        power[j] <- sum(D) / sum(hypo[, j])
      }
      FWER_Spending_c[l, b] <- mean(V)
      power_Spending_c[l, b] <- mean(power, na.rm = TRUE)


      ## ADDIS-Graph
      V <- rep(0, m) # Indicates, whether there was at least one type 1 error in a trial
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

  # Save results in data frame

  pi_A_vals <- rep(seq(0.1, 0.9, 0.1), times = length(batch_sizes))
  batch_size_vals <- rep(batch_sizes, each = 9)

  # Create data frames for power and FWER

  df_power <- data.frame(
    pi_A = rep(pi_A_vals, 3),
    batch_size = rep(batch_size_vals, 3),
    power = c(
      as.vector(power_Spending),
      as.vector(power_Spending_c),
      as.vector(power_Graph)
    ),
    procedure = rep(c("Spending", "Spending_c", "Graph"), each = length(pi_A_vals))
  )


  df_fwer <- data.frame(
    pi_A = rep(pi_A_vals, 3),
    batch_size = rep(batch_size_vals, 3),
    FWER = c(
      as.vector(FWER_Spending),
      as.vector(FWER_Spending_c),
      as.vector(FWER_Graph)
    ),
    procedure = rep(c("Spending", "Spending_c", "Graph"), each = length(pi_A_vals))
  )

  # Combine power and FWER data frames

  df_power$Metric <- "Power"
  df_fwer$Metric <- "FWER"
  colnames(df_fwer)[3] <- "value"
  colnames(df_power)[3] <- "value"

  df_combined <- rbind(df_power, df_fwer)

  # Save the data frame

  saveRDS(df_combined, file = filename)
}


### Function to generate the plots in the FWER case for different batch sizes

# Input:
# input_filename:     filename of the dataset with the results
# comparison_proc:    procedure that is compared with ADDIS-Graph, either "Spend" (ADDIS-Spending)
#                     or "Spend_c" (Closed ADDIS-Spending)

# Output: ggplot that compares the ADDIS-Graph with comparison_proc

plot_generator_FWER_batchsize <- function(input_filename, comparison_proc) {
  # Set alpha, as it is needed for horizontal line

  alpha <- 0.2

  # load the results

  results <- readRDS(input_filename)

  # ADDIS-Spending

  df_spending <- subset(results, procedure == comparison_proc)

  if (comparison_proc == "Spending") {
    p_spending <- ggplot(df_spending, aes(x = pi_A, y = value, linetype = factor(batch_size))) +
      geom_line(aes(color = Metric)) +
      geom_point(aes(color = Metric)) +
      geom_hline(yintercept = alpha, linetype = "dashed") +
      scale_linetype_manual(name = "Batch-size", values = c("1" = "solid", "5" = "longdash", 
                                                            "10" = "dashed", "20" = "dotted")) +
      scale_color_manual(guide = "none", values = c("FWER" = "#f84f4f", "Power" = "#f84f4f")) +
      scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
      xlab(expression(pi[A])) +
      ylab("FWER / Power") +
      ggtitle(expression("ADDIS-Spending"[local])) +
      theme_minimal() +
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
  } else {
    p_spending <- ggplot(df_spending, aes(x = pi_A, y = value, linetype = factor(batch_size))) +
      geom_line(aes(color = Metric)) +
      geom_point(aes(color = Metric)) +
      geom_hline(yintercept = alpha, linetype = "dashed") +
      scale_linetype_manual(name = "Batch-size", values = c("1" = "solid", "5" = "longdash", 
                                                            "10" = "dashed", "20" = "dotted")) +
      scale_color_manual(guide = "none", values = c("FWER" = "#f84f4f", "Power" = "#f84f4f")) +
      scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
      xlab(expression(pi[A])) +
      ylab("FWER / Power") +
      ggtitle(expression("Closed ADDIS-Spending"[local])) +
      theme_minimal() +
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
  }


  # ADDIS-Graph

  df_graph <- subset(results, procedure == "Graph")

  p_graph <- ggplot(df_graph, aes(x = pi_A, y = value, linetype = factor(batch_size))) +
    geom_line(aes(color = Metric)) +
    geom_point(aes(color = Metric, shape = Metric)) +
    geom_hline(yintercept = alpha, linetype = "dashed") +
    scale_linetype_manual(name = "Batch-size", values = c("1" = "solid", "5" = "longdash", 
                                                          "10" = "dashed", "20" = "dotted")) +
    scale_color_manual(guide = "none", values = c("FWER" = "cornflowerblue", "Power" = "cornflowerblue")) +
    scale_shape_manual(guide = "none", values = c("FWER" = 17, "Power" = 17)) +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    xlab(expression(pi[A])) +
    ylab("FWER / Power") +
    ggtitle(expression("ADDIS-Graph"[conf - u])) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

  # Combined

  combined_plot <- p_spending + p_graph & theme(legend.position = "bottom")
  combined_plot <- combined_plot + plot_layout(guides = "collect")

  return(combined_plot)
}


### Function to generate the data in the FWER case for time-varying batch-sizes

# Input:
# smallest_batches: smallest batch-sizes (some vector containing values from {1,2,5,10})
# direction:        specifies whether the batch-size is decreasing ("decr") or increasing ("incr")
# filename: Output filename

# Output: Saves simulation results (power and FWER) of the ADDIS-Spending and ADDIS-Graph
#         for the batch-sizes that vary over time in file "filename"

data_generator_FWER_smallest_batch <- function(smallest_batches, direction, filename) {
  ### Simulation parameters
  m <- 1000 # Number of Trials
  n <- 100
  mu_A <- 3 # Strength of the alternative
  mu_N <- -0.5 # conservativeness of null p-values
  corr <- 0.5 # correlation within batch
  # pi_A is defined in the loop below

  ### Initialise Hyperparameters
  alpha <- 0.2
  tau <- 0.8
  lambda <- 0.16
  gamma <- (6 / (pi^2 * (1:100)^2))

  ### Set seed to make the results reproducible
  set.seed(12345)

  ### Predefine vectors for FWER and power of the different procedures
  FWER_Spending <- matrix(0, 9, length(smallest_batches))
  power_Spending <- matrix(0, 9, length(smallest_batches))
  FWER_Graph <- matrix(0, 9, length(smallest_batches))
  power_Graph <- matrix(0, 9, length(smallest_batches))

  b <- 1 # Counter

  for (smallest_batch in smallest_batches) {
    numb_batch <- n / (10 * smallest_batch) # number each of the batch_sizes (s, 2s, 3s, 4s) is repeated,
    # where s is the smallest batch_size

    # Lags for the local dependence structure
    
    if(direction == "incr"){
      lags <- c(
        rep((0:(smallest_batch - 1)), numb_batch), rep((0:(2 * smallest_batch - 1)), numb_batch),
        rep((0:(3 * smallest_batch - 1)), numb_batch), rep((0:(4 * smallest_batch - 1)), numb_batch)
      )
    }else{
      lags <- c(
        rep((0:(4 * smallest_batch - 1)), numb_batch), rep((0:(3 * smallest_batch - 1)), numb_batch), 
        rep((0:(2 * smallest_batch - 1)), numb_batch), rep((0:(smallest_batch - 1)), numb_batch)
      )
    }
    

    # Calculate d_j for ADDIS-Graph
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
      # Fast way for p-values in batches.
      pi_A <- l / 10
      p <- matrix(, nrow = n, ncol = m)
      hypo <- matrix(, nrow = n, ncol = m)
      for (j in 1:m) {
        hypo[, j] <- rbinom(n, 1, pi_A)
        X <- rep(0, n)
        if(direction == "incr"){
          for (a in 1:4) {
            for (k in 1:numb_batch) {
              # parameters for correlation structure
              sigma <- matrix(corr, smallest_batch * a, smallest_batch * a) + diag((1 - corr), smallest_batch * a)
              mu <- rep(0, smallest_batch * a)

              place_X <- (smallest_batch * numb_batch * sum(1:(a - 1))) * (a > 1)
              # generation of test statistic
              X[(((k - 1) * smallest_batch * a + 1) + place_X):((k * smallest_batch * a) + place_X)] <- mvrnorm(1, mu, sigma)
            }
          }
        }else{
          for (a in 4:1) {
            for (k in 1:numb_batch) {
              # parameters for correlation structure
              sigma <- matrix(corr, smallest_batch * a, smallest_batch * a) + diag((1 - corr), smallest_batch * a)
              mu <- rep(0, smallest_batch * a)
              
              place_X <- (smallest_batch * numb_batch * sum(4:(a + 1))) * (a < 4)
              # generation of test statistic
              X[(((k - 1) * smallest_batch * a + 1) + place_X):((k * smallest_batch * a) + place_X)] <- mvrnorm(1, mu, sigma)
            }
          }
        }
        Z <- mu_N * (1 - hypo[, j]) + mu_A * hypo[, j] + X
        p[, j] <- pnorm(-Z)
      }

      ## ADDIS-Spending
      V <- rep(0, m) # Indicates, whether there was at least one type 1 error in a trial
      power <- rep(0, m) # Power within each trial

      for (j in 1:m) {
        alpha_ind <- ADDIS_Spending(alpha, gamma, tau, lambda, lags, p[, j], n)
        hypo_est <- alpha_ind >= p[, j]
        V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
        D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
        power[j] <- sum(D) / sum(hypo[, j])
      }
      FWER_Spending[l, b] <- mean(V)
      power_Spending[l, b] <- mean(power, na.rm = TRUE)


      ## ADDIS-Graph
      V <- rep(0, m) # Indicates, whether there was at least one type 1 error in a trial
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

  # Save results in data frame

  pi_A_vals <- rep(seq(0.1, 0.9, 0.1), times = length(smallest_batches))
  smallest_batches_vals <- rep(smallest_batches, each = 9)

  # Create data frames for power and FWER

  df_power <- data.frame(
    pi_A = rep(pi_A_vals, 2),
    smallest_batches = rep(smallest_batches_vals, 2),
    power = c(
      as.vector(power_Spending),
      as.vector(power_Graph)
    ),
    procedure = rep(c("Spending", "Graph"), each = length(pi_A_vals))
  )


  df_fwer <- data.frame(
    pi_A = rep(pi_A_vals, 2),
    smallest_batches = rep(smallest_batches_vals, 2),
    FWER = c(
      as.vector(FWER_Spending),
      as.vector(FWER_Graph)
    ),
    procedure = rep(c("Spending", "Graph"), each = length(pi_A_vals))
  )

  # Combine power and FWER data frames

  df_power$Metric <- "Power"
  df_fwer$Metric <- "FWER"
  colnames(df_fwer)[3] <- "value"
  colnames(df_power)[3] <- "value"

  df_combined <- rbind(df_power, df_fwer)

  # Save the data frame

  saveRDS(df_combined, file = filename)
}


### Function to generate the plots in the FWER case for time-varying batch sizes

# Input:
# input_filename:     filename of the dataset with the results

# Output: ggplot that compares the ADDIS-Graph with ADDIS-Spending using the results in input_filename

plot_generator_FWER_smallest_batch <- function(input_filename) {
  # Set alpha, as it is needed for horizontal line

  alpha <- 0.2

  # load the results

  results <- readRDS(input_filename)

  # ADDIS-Spending

  df_spending <- subset(results, procedure == "Spending")

  p_spending <- ggplot(df_spending, aes(x = pi_A, y = value, linetype = factor(smallest_batches))) +
    geom_line(aes(color = Metric)) +
    geom_point(aes(color = Metric)) +
    geom_hline(yintercept = alpha, linetype = "dashed") +
    scale_linetype_manual(name = "Smallest batch-size", values = c("1" = "solid", "2" = "longdash", 
                                                                   "5" = "dashed", "10" = "dotted")) +
    scale_color_manual(guide = "none", values = c("FWER" = "#f84f4f", "Power" = "#f84f4f")) +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    xlab(expression(pi[A])) +
    ylab("FWER / Power") +
    ggtitle(expression("ADDIS-Spending"[local])) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))


  # ADDIS-Graph

  df_graph <- subset(results, procedure == "Graph")

  p_graph <- ggplot(df_graph, aes(x = pi_A, y = value, linetype = factor(smallest_batches))) +
    geom_line(aes(color = Metric)) +
    geom_point(aes(color = Metric, shape = Metric)) +
    geom_hline(yintercept = alpha, linetype = "dashed") +
    scale_linetype_manual(name = "Smallest batch-size", values = c("1" = "solid", "2" = "longdash", 
                                                                   "5" = "dashed", "10" = "dotted")) +
    scale_color_manual(guide = "none", values = c("FWER" = "cornflowerblue", "Power" = "cornflowerblue")) +
    scale_shape_manual(guide = "none", values = c("FWER" = 17, "Power" = 17)) +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    xlab(expression(pi[A])) +
    ylab("FWER / Power") +
    ggtitle(expression("ADDIS-Graph"[conf - u])) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

  # Combined

  combined_plot <- p_spending + p_graph & theme(legend.position = "bottom")
  combined_plot <- combined_plot + plot_layout(guides = "collect")

  return(combined_plot)
}



### Function to generate the data in the FDR case for different test durations

# Input:
# gamma:    initial weights (nonnegative n-dim. vector)
# filename: Output filename

# Output: Saves simulation results (power and FDR) of the ADDIS* and FDR-ADDIS-Graph for the different
#         test durations (0, 1, 2, 5) of the asynchrony structure in file "filename"

data_generator_FDR_async <- function(gamma, filename) {
  ### Simulation parameters

  m <- 1000 # Number of Trials
  n <- 100 # Number of Hypotheses per Trial
  mu_A <- 3 # Strength of the alternative
  mu_N <- -0.5 # Conservativeness of null p-values (<0 for conservative null p-values)
  # pi_A is defined in the loop below

  ### Initialise Hyperparameters
  alpha <- 0.05
  tau <- 0.5
  lambda <- 0.25

  # Set stopping times
  e_values <- c(0, 1, 2, 5)


  ### Predefine vectors for FWER and power of the different procedures
  FDR_ADDIS <- matrix(0, 9, length(e_values))
  power_ADDIS <- matrix(0, 9, length(e_values))
  FDR_Graph <- matrix(0, 9, length(e_values))
  power_Graph <- matrix(0, 9, length(e_values))

  ### Set seed to make the results reproducible
  set.seed(12345)

  b <- 1 # Counter
  for (e_value in e_values) {
    e <- seq(1, n) + rep(e_value, n)

    # Adjusted weigths
    w_adj <- matrix(0, n, n)

    for (k in 1:n) {
      if (e[k] >= n) {
        break
      } else if (e[k] == k) {
        w_adj[k, (k + 1):n] <- gamma[1:(n - k)]
      } else {
        gamma_exh <- gamma[(e[k] - k + 1):(n - k)] / (1 - sum(gamma[1:(e[k] - k)]))
        w_adj[k, (e[k] + 1):n] <- gamma_exh
      }
    }


    ### Generate p-values and compute FDR and power for the desired procedures
    for (l in 1:9) {
      pi_A <- l / 10
      p <- matrix(, nrow = n, ncol = m)
      hypo <- matrix(, nrow = n, ncol = m)
      for (j in 1:m) {
        hypo[, j] <- rbinom(n, 1, pi_A)
        X <- rnorm(n)
        Z <- mu_N * (hypo[, j] - 1) * (-1) + mu_A * hypo[, j] + X
        p[, j] <- pnorm(-Z)
      }

      ## ADDIS*
      power <- rep(0, m) # Power within each trial
      fdr <- rep(0, m) # FDR within each trial

      for (j in 1:m) {
        alpha_ind <- ADDIS_star(alpha, gamma, tau, lambda, e, p[, j], n)
        hypo_est <- alpha_ind >= p[, j]
        D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
        power[j] <- sum(D) / sum(hypo[, j])
        fdr[j] <- sum((hypo[, j] == 0 & hypo_est == 1)) / max(1, sum(hypo_est))
      }
      FDR_ADDIS[l, b] <- mean(fdr)
      power_ADDIS[l, b] <- mean(power)


      ## ADDIS-Graph FDR
      power <- rep(0, m) # Power within each trial
      fdr <- rep(0, m) # FDR within each trial

      for (j in 1:m) {
        alpha_ind <- ADDIS_Graph_fdr(alpha, gamma, w_adj, w_adj, tau, lambda, e, p[, j], n)
        hypo_est <- alpha_ind >= p[, j]
        D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
        power[j] <- sum(D) / sum(hypo[, j])
        fdr[j] <- sum((hypo[, j] == 0 & hypo_est == 1)) / max(1, sum(hypo_est))
      }
      FDR_Graph[l, b] <- mean(fdr)
      power_Graph[l, b] <- mean(power)
    }
    b <- b + 1
  }

  # Save results in data frame

  pi_A_vals <- rep(seq(0.1, 0.9, 0.1), times = length(e_values))
  duration_vals <- rep(e_values, each = 9)

  # Create data frames for power and FDR

  df_power <- data.frame(
    pi_A = rep(pi_A_vals, 2),
    test_duration = rep(duration_vals, 2),
    power = c(
      as.vector(power_ADDIS),
      as.vector(power_Graph)
    ),
    procedure = rep(c("ADDIS", "Graph"), each = length(pi_A_vals))
  )


  df_fdr <- data.frame(
    pi_A = rep(pi_A_vals, 2),
    test_duration = rep(duration_vals, 2),
    power = c(
      as.vector(FDR_ADDIS),
      as.vector(FDR_Graph)
    ),
    procedure = rep(c("ADDIS", "Graph"), each = length(pi_A_vals))
  )

  # Combine power and FDR data frames

  df_power$Metric <- "Power"
  df_fdr$Metric <- "FDR"
  colnames(df_fdr)[3] <- "value"
  colnames(df_power)[3] <- "value"

  df_combined <- rbind(df_power, df_fdr)

  # Save the data frame

  saveRDS(df_combined, file = filename)
}


### Function to generate the plots in the FDR case for different test durations

# Input:
# input_filename:     filename of the dataset with the results

# Output: ggplot that compares the FDR-ADDIS-Graph with ADDIS*

plot_generator_FDR_async <- function(input_filename) {
  # Set alpha, as it is needed for horizontal line

  alpha <- 0.05

  # load the results

  results <- readRDS(input_filename)

  # ADDIS-Spending

  df_spending <- subset(results, procedure == "ADDIS")

  p_spending <- ggplot(df_spending, aes(x = pi_A, y = value, linetype = factor(test_duration))) +
    geom_line(aes(color = Metric)) +
    geom_point(aes(color = Metric)) +
    geom_hline(yintercept = alpha, linetype = "dashed") +
    scale_linetype_manual(name = "Test duration", values = c("0" = "solid", "1" = "longdash", 
                                                             "2" = "dashed", "5" = "dotted")) +
    scale_color_manual(guide = "none", values = c("FDR" = "#f84f4f", "Power" = "#f84f4f")) +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    xlab(expression(pi[A])) +
    ylab("FDR / Power") +
    ggtitle(expression("ADDIS*"[async])) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))



  # ADDIS-Graph

  df_graph <- subset(results, procedure == "Graph")

  p_graph <- ggplot(df_graph, aes(x = pi_A, y = value, linetype = factor(test_duration))) +
    geom_line(aes(color = Metric)) +
    geom_point(aes(color = Metric, shape = Metric)) +
    geom_hline(yintercept = alpha, linetype = "dashed") +
    scale_linetype_manual(name = "Test duration", values = c("0" = "solid", "1" = "longdash", 
                                                             "2" = "dashed", "5" = "dotted")) +
    scale_color_manual(guide = "none", values = c("FDR" = "cornflowerblue", "Power" = "cornflowerblue")) +
    scale_shape_manual(guide = "none", values = c("FDR" = 17, "Power" = 17)) +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    xlab(expression(pi[A])) +
    ylab("FDR / Power") +
    ggtitle(expression("FDR-ADDIS-Graph"[async])) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

  # Combined

  combined_plot <- p_spending + p_graph & theme(legend.position = "bottom")
  combined_plot <- combined_plot + plot_layout(guides = "collect")

  return(combined_plot)
}


### Function to generate the data in the FWER case with incorporation of the correlation structure
### for 4 different batch-sizes / correlation strengths / conservativeness of null p-values

# Input:
# mu_Ns:    conservativeness of null p-values (4-dimensional vector with values <=0)
# rhos:     correlation within batch (4-dimensional vector with values in [0,1))
# batch_sizes:    batch-sizes (4-dimensional vector with integers)
# var_ind:    indicator where mu_Ns, rhos or batch-sizes should vary (choose one of ("mu_N", "rho", "batch"))
# filename: Output filename

# Output: Saves simulation results (power and FWER) of the Adaptive_Graph_{corr} and
#         ADDIS-Graph_{conf} for different values of var_ind in file "filename"

data_generator_FWER_corr <- function(mu_Ns = c(), rhos = c(), batch_sizes = c(), var_ind, filename) {
  ### Simulation parameters
  m <- 1000 # Number of Trials
  n <- 100 # Number of hypotheses
  mu_A <- 3 # Strength of the alternative
  mu_N <- 0 # Conservativeness of null p-values
  rho <- 0.5 # Correlation within batch
  batch_size <- 10 # batch-size
  # pi_A is defined in the loop below

  ### Initialise Hyperparameters
  alpha <- 0.2
  gamma <- (6 / (pi^2 * (1:100)^2))

  if (var_ind == "mu_N") {
    tau <- 0.8
    lambda <- 0.16
  } else {
    tau <- 1
    lambda <- 0.2
  }

  lambda_corr <- 0.2

  ### Set seed to make the results reproducible
  set.seed(12345)

  ### Predefine vectors for FWER and power of the different procedures
  FWER_corr <- matrix(0, 9, 4)
  power_corr <- matrix(0, 9, 4)
  FWER_Graph <- matrix(0, 9, 4)
  power_Graph <- matrix(0, 9, 4)

  for (b in 1:4) {
    if (var_ind == "mu_N") {
      mu_N <- mu_Ns[b]
    } else if (var_ind == "rho") {
      rho <- rhos[b]
    } else {
      batch_size <- batch_sizes[b]
    }

    ### Parameters for a local dependence structure given by batches
    batch_number <- ceiling(n / batch_size)
    lags <- rep(seq(0, (batch_size - 1), 1), batch_number)
    lags <- lags[1:n]
    sigma <- matrix(rho, batch_size, batch_size) + diag((1 - rho), batch_size)
    mu <- rep(0, batch_size)
    corr <- as.matrix(bdiag(replicate(batch_number, sigma, simplify = FALSE)))

    # Calculate d_j for ADDIS-Graph
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

    # Calculate local dependence adjusted weights
    w_lda <- matrix(0, n, n)
    for (k in 1:n) {
      if (d_j[k] > n) {
        break
      } else if (d_j[k] == k + 1) {
        w_lda[k, d_j[k]:n] <- gamma[1:(n - k)]
      } else {
        gamma_exh <- gamma[(d_j[k] - k):(n - k)] / (1 - sum(gamma[1:(d_j[k] - k - 1)]))
        w_lda[k, d_j[k]:n] <- gamma_exh
      }
    }

    ### Generate p-values and compute FWER and power for the desired procedures
    for (l in 1:9) {
      ## Fast way for p-values in batches.
      pi_A <- l / 10
      p <- matrix(, nrow = n, ncol = m)
      Z <- matrix(, nrow = n, ncol = m)
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

      ## Adaptive_graph_{corr}
      V <- rep(0, m) # Indicates, whether there was atleast one type 1 error in a trial
      power <- rep(0, m) # Power within each trial

      for (j in 1:m) {
        alpha_ind <- ADDIS_Graph_corr(alpha, corr, gamma, w_lda, lambda_corr, lags, p[, j], n)$alpha_ind
        hypo_est <- alpha_ind >= p[, j]
        V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
        D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
        power[j] <- sum(D) / sum(hypo[, j])
      }
      FWER_corr[l, b] <- mean(V)
      power_corr[l, b] <- mean(power, na.rm = TRUE)


      ## ADDIS-Graph_{local-u}
      V <- rep(0, m) # Indicates, whether there was atleast one type 1 error in a trial
      power <- rep(0, m) # Power within each trial

      for (j in 1:m) {
        alpha_ind <- ADDIS_Graph(alpha, gamma, w_lda, tau, lambda, lags, p[, j], n)
        hypo_est <- alpha_ind >= p[, j]
        V[j] <- max((hypo[, j] == 0 & hypo_est == 1))
        D <- (hypo[, j] == 1 & hypo[, j] == hypo_est)
        power[j] <- sum(D) / sum(hypo[, j])
      }
      FWER_Graph[l, b] <- mean(V)
      power_Graph[l, b] <- mean(power, na.rm = TRUE)
    }
  }

  # Save results in data frame

  pi_A_vals <- rep(seq(0.1, 0.9, 0.1), times = 4)

  if (var_ind == "mu_N") {
    var_vals <- rep(mu_Ns, each = 9)
  } else if (var_ind == "rho") {
    var_vals <- rep(rhos, each = 9)
  } else {
    var_vals <- rep(batch_sizes, each = 9)
  }

  # Create data frames for power and FWER

  df_power <- data.frame(
    pi_A = rep(pi_A_vals, 2),
    var = rep(var_vals, 2),
    power = c(
      as.vector(power_corr),
      as.vector(power_Graph)
    ),
    procedure = rep(c("Corr", "Graph"), each = length(pi_A_vals))
  )


  df_fwer <- data.frame(
    pi_A = rep(pi_A_vals, 2),
    var = rep(var_vals, 2),
    power = c(
      as.vector(FWER_corr),
      as.vector(FWER_Graph)
    ),
    procedure = rep(c("Corr", "Graph"), each = length(pi_A_vals))
  )

  # Combine power and FWER data frames

  df_power$Metric <- "Power"
  df_fwer$Metric <- "FWER"
  colnames(df_fwer)[3] <- "value"
  colnames(df_power)[3] <- "value"

  df_combined <- rbind(df_power, df_fwer)

  # Save the data frame

  saveRDS(df_combined, file = filename)
}


### Function to generate the plots in the FWER case with incorporation of the correlation structure
### for 4 different batch-sizes / correlation strengths / conservativeness of null p-values

# Input:
# input_filename:   filename of the dataset with the results
# var_ind:          the variable that is varied, choose from ("mu_N", "rho", "batch")

# Output: ggplot that compares the Adaptive-Graph_{corr} with ADDIS-Graph_{conf}/Adaptive-Graph_{conf}

plot_generator_FWER_corr <- function(input_filename, var_ind) {
  # Set alpha, as it is needed for horizontal line

  alpha <- 0.2

  # load the results

  results <- readRDS(input_filename)

  # With correlation

  df_corr <- subset(results, procedure == "Corr")

  p_corr <- ggplot(df_corr, aes(x = pi_A, y = value, linetype = factor(var))) +
    geom_line(aes(color = Metric)) +
    geom_point(aes(color = Metric)) +
    geom_hline(yintercept = alpha, linetype = "dashed") +
    scale_linetype_manual(guide = "none") +
    scale_color_manual(guide = "none", values = c("FWER" = "#ff5c33", "Power" = "#ff5c33")) +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    xlab(expression(pi[A])) +
    ylab("FWER / Power") +
    ggtitle(expression("ADDIS-Spending"[local])) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

  if (var_ind == "batch") {
    p_corr <- p_corr + scale_linetype_manual(name = "Batch-size", values = c("1" = "solid", "5" = "longdash", 
                                                                             "10" = "dashed", "20" = "dotted")) +
      ggtitle(expression("Adaptive-Graph"[corr]))
  } else if (var_ind == "mu_N") {
    p_corr <- p_corr + scale_linetype_manual(name = expression(paste(mu[N], "  ")), values = c("0" = "solid", "-0.5" = "longdash", 
                                                                                               "-1" = "dashed", "-2" = "dotted")) +
      ggtitle(expression("Adaptive-Graph"[corr]))
  } else {
    p_corr <- p_corr + scale_linetype_manual(name = expression(paste(rho, "  ")), values = c("0.3" = "solid", "0.5" = "longdash", 
                                                                                             "0.7" = "dashed", "0.9" = "dotted")) +
      ggtitle(expression("Adaptive-Graph"[corr]))
  }

  # ADDIS-Graph

  df_graph <- subset(results, procedure == "Graph")

  p_graph <- ggplot(df_graph, aes(x = pi_A, y = value, linetype = factor(var))) +
    geom_line(aes(color = Metric)) +
    geom_point(aes(color = Metric, shape = Metric)) +
    geom_hline(yintercept = alpha, linetype = "dashed") +
    scale_linetype_manual(guide = "none") +
    scale_color_manual(guide = "none", values = c("FWER" = "cornflowerblue", "Power" = "cornflowerblue")) +
    scale_shape_manual(guide = "none", values = c("FWER" = 17, "Power" = 17)) +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    xlab(expression(pi[A])) +
    ylab("FWER / Power") +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

  if (var_ind == "batch") {
    p_graph <- p_graph + scale_linetype_manual(name = "Batch-size", 
                                               values = c("1" = "solid", "5" = "longdash", 
                                                          "10" = "dashed", "20" = "dotted")) +
      ggtitle(expression("Adaptive-Graph"[conf]))
  } else if (var_ind == "mu_N") {
    p_graph <- p_graph + scale_linetype_manual(name = expression(paste(mu[N], "  ")), 
                                               values = c("0" = "solid", "-0.5" = "longdash", 
                                                          "-1" = "dashed", "-2" = "dotted")) +
      ggtitle(expression("ADDIS-Graph"[conf]))
  } else {
    p_graph <- p_graph + scale_linetype_manual(name = expression(paste(rho, "  ")), 
                                               values = c("0.3" = "solid", "0.5" = "longdash", 
                                                          "0.7" = "dashed", "0.9" = "dotted")) +
      ggtitle(expression("Adaptive-Graph"[conf]))
  }

  # Combined

  combined_plot <- p_corr + p_graph & theme(legend.position = "bottom")
  combined_plot <- combined_plot + plot_layout(guides = "collect")

  return(combined_plot)
}
