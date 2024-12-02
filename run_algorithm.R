### Import and set directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  path <- getwd()
  print(getwd())
} else {
  message("The rstudioapi package is not available.")
}

source("required_packages.R")
source("Functions/support_functions.R")
source("Functions/initialization_functions.R")
source("Functions/E-step.R")
source("Functions/M-step.R")
source("Functions/log_likelihood.R")
source("Functions/EM-algorithm.R")
source("Functions/BIC_calculation.R")
source("Functions/table_BIC_comparison.R")
source("Functions/Predictions.R")
 
# INPUT:

# - data: complete dataset with response and all the covariates
# - C: fixed number of latent clusters
# - cont_var (U): data containing continuous covariates
# - cat_var (V): list containing independent categorical covariates
# - dich_dep_var (D): data containing dichotomous dependent covariates

# init can be "kmeans" or "random"

# mstep_isingfit = F: uses IsingSampler::EstimateIsing function to estimate the Ising Model
# Estimates the Ising model by maximizing the pseudolikelihood

# mstep_isingfit = T: uses IsingFit::IsingFit function to estimate the Ising Model
# Estimates the Ising model using the eLasso method, which combines l1-regularized
# logistic regression with model selection based on the Extended Bayesian Information Criterion (EBIC)

# max_it and tol: max_it and tolearnce for stopping the EM-algorithm

# Load simulated dataset with C = 3
data = read.csv( "Example Datasets/data1.csv" )

# Set number of clusters
C = 3

# Run the algorithm with N completely random initializations
N <- 25
log_ll <- vector(mode="numeric", length = N)
fitting_rec <- list()

set.seed(7)
for(i in 1:N){
  tryCatch({
    fitting <- fitMLCWMd(data = data, formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
                         C = C, 
                         cont_var =  c("x1", "x2"), 
                         cat_var = c("V1"), 
                         dich_dep_var = c("D1", "D2", "D3"), 
                         init = "random", 
                         mstep_isingfit = TRUE, 
                         fixed_intercept = TRUE,
                         max_it = 30, 
                         tol = 1e-6)
    log_ll[i] <- fitting$log_l[length(fitting$log_l)]
    fitting_rec[[i]] <- fitting
  }, error = function(err) {
    cat(paste0(
      "==============================\n",
      "Bad initialization detected\n",
      "------------------------------\n",
      "Index: ", i, "\n",
      "Error Message: ", conditionMessage(err), "\n",
      "Attempting next initialization...\n",
      "==============================\n"
    ))
  }
  )
}

# Select the result with the highest log-likelihood
log_ll <- if_else(log_ll == 0, -Inf, log_ll)
fitting <- fitting_rec[[which.max(log_ll)]]
fitting$bic

# View Log-Likelihood trend
ggplot(data = data.frame(log_l = fitting$log_l[2:length(fitting$log_l)], iteration = 1:(length(fitting$log_l)-1)), 
       aes(x = iteration, y = log_l)) + 
  geom_line(color = 'orange') +
  geom_point(color = 'red', size = 3, pch = 20) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 100, by = 1)) +
  labs(x = "Iteration", y = "Log-likelihood")


# Get Predictions for a test data
# Load Test data
test_data = read.csv( "Example Datasets/test_data1.csv" )
y_pred <- predict_MLCWMd(fitting, test_data, C, new_level = FALSE)

table_Pred(y_pred, test_data, fitting)

# Run the algorithm with N k-means initializations
set.seed(7)
N <- 1
log_ll <- vector(mode="numeric", length = N)
fitting_rec <- list()

for(i in 1:N){
  fitting <- tryCatch({
    fitting <- fitMLCWMd(data = data, formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
                         C = C, 
                         cont_var =  c("x1", "x2"), 
                         cat_var = c("V1"), 
                         dich_dep_var = c("D1", "D2", "D3"), 
                         init = "kmeans", 
                         mstep_isingfit = TRUE, 
                         fixed_intercept = TRUE,
                         max_it = 30, 
                         tol = 1e-6)
    log_ll[i] <- fitting$log_l[length(fitting$log_l)]
    fitting_rec[[i]] <- fitting
  }, error = function(err) {
    log_ll[i] <- -Inf
    cat(paste0(
      "==============================\n",
      "Bad initialization detected\n",
      "------------------------------\n",
      "Index: ", i, "\n",
      "Error Message: ", conditionMessage(err), "\n",
      "Attempting next initialization...\n",
      "==============================\n"
    ))
  }
  )
}

# Select the result with the highest log-likelihood
log_ll <- if_else(log_ll == 0, -Inf, log_ll)
fitting <- fitting_rec[[which.max(log_ll)]]
fitting$bic


# View Log-Likelihood trend
ggplot(data = data.frame(log_l = fitting$log_l[2:length(fitting$log_l)], iteration = 1:(length(fitting$log_l)-1)), 
       aes(x = iteration, y = log_l)) + 
  geom_line(color = 'orange') +
  geom_point(color = 'red', size = 3, pch = 20) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 100, by = 1)) +
  labs(x = "Iteration", y = "Log-likelihood")

# Get Predictions for a test data
# Load Test data
test_data = read.csv( "Example Datasets/test_data1.csv" )
y_pred <- predict_MLCWMd(fitting, test_data, C, new_level = FALSE)

table_Pred(y_pred, test_data, fitting)

#--------------------Complete Procedure Described in the paper------------------

# PROCEDURE:
# - Set the assumed range of values for C.
# - For each C, run the algorithm 30 times with random initialization 
#   and select the result with the highest log-likelihood.
#   If all runs fail, switch to k-means initialization.
# - Calculate the BIC for the best iteration among all C values.
# - Finally, select the number of clusters with the lowest BIC.

# NOTE: Many iterations may result in errors. For instance, starting with a random 
# initialization that produces an extremely small cluster can lead to failure in 
# the mixed-effects regression estimation. If all random initializations fail, 
# we switch to k-means initialization.

set.seed(7)
possible_C_values <- c(1, 2, 3, 4)
saved_fitting <- list()
BIC_vector <- vector(mode = "numeric", length = length(possible_C_values))

for (k in 1:length(possible_C_values)) {
  fitting_rec <- list()
  N <- 25
  log_ll <- vector(mode = "numeric", length = N)
  success <- FALSE  # Track whether a successful fit occurs
  
  cat(paste("==============================\n",
            "Running with C = ", possible_C_values[k], "\n",
            "==============================\n"))
  
  for (i in 1:N) {
    tryCatch({
      fitting <- fitMLCWMd(data = data, 
                           formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
                           C = possible_C_values[k], 
                           cont_var = c("x1", "x2"), 
                           cat_var = c("V1"), 
                           dich_dep_var = c("D1", "D2", "D3"), 
                           init = "random", 
                           mstep_isingfit = TRUE, 
                           fixed_intercept = TRUE,
                           max_it = 30, 
                           tol = 1e-6)
      log_ll[i] <- fitting$log_l[length(fitting$log_l)]
      fitting_rec[[i]] <- fitting
      success <- TRUE  # Mark that a successful fit occurred
    }, error = function(err) {
    cat(paste0(
        "==============================\n",
        "Bad initialization detected\n",
        "------------------------------\n",
        "Initialization: random", "\n",
        "Number of clusters: ", k, "\n",
        "Index: ", i, "\n",
        "Error Message: ", conditionMessage(err), "\n",
        "Attempting next initialization...\n",
        "==============================\n"
      ))
    })
  }
  
  # If no successful runs, try a different initialization method
  if (!success) {
    cat(paste("All runs failed for C =", possible_C_values[k], "Retrying with a different initialization...\n"))
    for (i in 1:5) {
      tryCatch({
        fitting <- fitMLCWMd(data = data, 
                             formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
                             C = possible_C_values[k], 
                             cont_var = c("x1", "x2"), 
                             cat_var = c("V1"), 
                             dich_dep_var = c("D1", "D2", "D3"), 
                             init = "kmeans",  
                             mstep_isingfit = TRUE, 
                             fixed_intercept = TRUE,
                             max_it = 30, 
                             tol = 1e-6)
        log_ll[i] <- fitting$log_l[length(fitting$log_l)]
        fitting_rec[[i]] <- fitting
        success <- TRUE
      }, error = function(err) {
        cat(paste0(
          "==============================\n",
          "Bad initialization detected\n",
          "------------------------------\n",
          "Initialization: kmeans", "\n",
          "Number of clusters: ", k, "\n",
          "Index: ", i, "\n",
          "Error Message: ", conditionMessage(err), "\n",
          "Attempting next initialization...\n",
          "==============================\n"
        ))
      })
    }
  }
  
  # Handle log_ll values and select the best fitting
  log_ll <- if_else(log_ll == 0, -Inf, log_ll)
  if (any(is.finite(log_ll))) {
    best_fit_index <- which.max(log_ll)
    fitting <- fitting_rec[[best_fit_index]]
    BIC_vector[k] <- fitting$bic
    saved_fitting[[k]] <- fitting
  } else {
    cat(paste("No successful fits for C =", possible_C_values[k], "\n"))
    BIC_vector[k] <- NA
    saved_fitting[[k]] <- NULL
  }
}


# See the best number of clusters identified by the BIC
table_BIC(BIC_vector, possible_C_values)

# Best Fit
best_fit <- saved_fitting[[3]]
best_fit$bic

# View Log-Likelihood trend
ggplot(data = data.frame(log_l = fitting$log_l[2:length(best_fit$log_l)], 
                         iteration = 1:(length(best_fit$log_l)-1)), 
       aes(x = iteration, y = log_l)) + 
  geom_line(color = 'orange') +
  geom_point(color = 'red', size = 3, pch = 20) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 100, by = 1)) +
  labs(x = "Iteration", y = "Log-likelihood")

# Get Predictions for a test data
# Load Test data
test_data = read.csv( "Example Datasets/test_data1.csv" )
y_pred <- predict_MLCWMd(fitting, test_data, C, new_level = FALSE)

table_Pred(y_pred, test_data, fitting)
