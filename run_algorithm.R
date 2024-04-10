source("required_packages.R")
source("Functions/support_functions.R")
source("Functions/initialization_functions.R")
source("Functions/E-step.R")
source("Functions/M-step.R")
source("Functions/log_likelihood.R")
source("Functions/EM-algorithm.R")
source("Functions/BIC_calculation.R")
source("Functions/table_BIC_comparison.R")
source("Example Datasets/example_dataset_1.R")
#source("Example Datasets/example_dataset_2.R")
 

# INPUT:

# - data: complete dataset with response and all the covariates
# - C: fixed number of latent clusters
# - U: data containing continuous covariates
# - V: list containing independent categorical covariates
# - D: data containing dichotomous dependent covariates

# flag_init = T: uses completely random initialization
# flag_init = F: uses k-means initialization

# flag_mstep = T: uses IsingSampler::EstimateIsing function to estimate the Ising Model
# Estimates the Ising model by maximizing the pseudolikelihood

# flag_mstep = F: uses IsingFit::IsingFit function to estimate the Ising Model
# Estimates the Ising model using the eLasso method, which combines l1-regularized
# logistic regression with model selection based on the Extended Bayesian Information Criterion (EBIC)

# Set number of clusters

C = 3

# Run the algorithm with N completely random initializations
N <- 30
log_ll <- vector(mode="numeric", length = N)
fitting_rec <- list()

for(i in 1:N){
  tryCatch({
    fitting <- Algo_full(data = data, C = C, U = U, V = V, D = D, flag_init = T, flag_mstep = F)
    log_ll[i] <- fitting$log_l[length(fitting$log_l)]
    fitting_rec[[i]] <- fitting
  }, error = function(err) {
    cat(paste("Error occurred at index", i, "Error message:", conditionMessage(err), 
              "\n", "Next Initialization...", "\n"))
  }
  )
}

# Select the result with the highest log-likelihood
log_ll <- if_else(log_ll == 0, -Inf, log_ll)
fitting <- fitting_rec[[which.max(log_ll)]]

# Run the algorithm with N k-means initializations
N <- 10
log_ll <- vector(mode="numeric", length = N)
fitting_rec <- list()

for(i in 1:N){
  fitting <- tryCatch({
    fitting <- Algo_full(data = data, C = C, U = U, V = V, D = D, flag_init = F, flag_mstep = F)
    log_ll[i] <- fitting$log_l[length(fitting$log_l)]
    fitting_rec[[i]] <- fitting
  }, error = function(err) {
    log_ll[i] <- -Inf
    cat(paste("Error occurred at index", i, "Error message:", conditionMessage(err), 
              "\n", "Next Initialization...", "\n"))
  }
  )
}

# Select the result with the highest log-likelihood
log_ll <- if_else(log_ll == 0, -Inf, log_ll)
fitting <- fitting_rec[[which.max(log_ll)]]


# View Log-Likelihood trend
ggplot(data = data.frame(log_l = fitting$log_l[2:length(fitting$log_l)], iteration = 1:(length(fitting$log_l)-1)), 
       aes(x = iteration, y = log_l)) + 
  geom_line(color = 'orange') +
  geom_point(color = 'red', size = 3, pch = 20) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 100, by = 1)) +
  labs(x = "Iteration", y = "Log-likelihood")

#--------------------Complete Procedure Described in the paper------------------

# - Set the assumed range of values for C

# - For each c, run the algorithm 30 times with random initialization 
#   and select the one with the highest loglik

# - Calculate the BIC for all the best iteration among all c values

# - At the end, select the number of cluster with the lowest BIC

# NOTE: Many of the iterations will result in errors. For instance, if we start, 
# by taking a random initialization, with an extremely small cluster, it may lead 
# to failure in the mixed effects regression estimation.

possible_C_values <- c(2,3,4)
saved_fitting <- list()
BIC_vector <- vector(mode="numeric", length = length(possible_C_values))
for(k in 1:length(possible_C_values)){
 fitting_rec <- list()
 N <- 30
 log_ll <- vector(mode="numeric", length = N)
 for(i in 1:N){
  tryCatch({
    fitting <- Algo_full(data = data, C = possible_C_values[k], U = U, V = V, D = D, 
                         flag_init = F, flag_mstep = F)
    log_ll[i] <- fitting$log_l[length(fitting$log_l)]
    fitting_rec[[i]] <- fitting
  }, error = function(err) {
    cat(paste("Error occurred at index", i,  "for number of clusters C = ", k+1, "|","Error message:", conditionMessage(err), 
              "\n", "Next Initialization...", "\n \n"))
  }
  )
 }
 log_ll <- if_else(log_ll == 0, -Inf, log_ll)
 fitting <- fitting_rec[[which.max(log_ll)]]

 BIC_vector[k] <- Calculate_BIC(U, V, D, 
                                log_l = fitting$log_l[length(fitting$log_l)], possible_C_values[k], 
                                fixed_intercept = FALSE)
 saved_fitting[[k]] <- fitting 
 cat(paste("Starting running with C =", possible_C_values[k+1], "\n \n"))
}

# See the best number of clusters identified by the BIC
table_BIC(vec, possible_C_values)
