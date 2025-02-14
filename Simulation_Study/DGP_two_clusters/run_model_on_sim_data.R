### Import and set directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  path <- getwd()
  print(getwd())
} else {
  message("The rstudioapi package is not available.")
}

source( paste0(dirname(dirname(path)), "/required_packages.R") )
source( paste0(dirname(dirname(path)), "/Functions/EM-algorithm.R") )
source( paste0(dirname(dirname(path)), "/Functions/support_functions.R") )
source( paste0(dirname(dirname(path)), "/Functions/initialization_functions.R") )
source( paste0(dirname(dirname(path)), "/Functions/E-step.R") )
source( paste0(dirname(dirname(path)), "/Functions/M-step.R") )
source( paste0(dirname(dirname(path)), "/Functions/log_likelihood.R") )
source( paste0(dirname(dirname(path)), "/Functions/BIC_calculation.R") )


###### FIT THE MODEL #########

fit_list_full = vector("list",3)
fit_list_partial = vector("list",3)

N <- 30
R <- 100
dir.create("fitted_models")

for (i in 1:R){
  
  set.seed(i*7)
  
  data = readRDS( paste0("sim_data/mydata", i, ".rds") )
  
  for (c in 1:3){
    
    ### Initialize lists to store BICs and models
    bic_full_list <- numeric(N)
    bic_partial_list <- numeric(N)
    model_full_list <- vector("list", N)
    model_partial_list <- vector("list", N)
    
    ### Run each algorithm N times and compute BIC for full model
    for (n in 1:N) {
      fit_c_full = tryCatch(
        fitMLCWMd(data = data, formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
                  C = c, 
                  cont_var =  c("x1", "x2"), 
                  cat_var = c("V1"), 
                  dich_dep_var = c("D1", "D2", "D3"), 
                  init = "kmeans", 
                  mstep_isingfit = TRUE, ###
                  fixed_intercept = TRUE,
                  max_it = 30, 
                  tol = 1e-6),
        error = function(e) {
          list(model = NA)
        }
      )
      
      # Calculate BIC if model is not NA
      if (is.null(fit_c_full$model)) {
        bic_full_list[n] <- fit_c_full$bic
      } else {
        bic_full_list[n] <- Inf  # Set to infinity if model fails
      }
      
      model_full_list[[n]] <- fit_c_full
    }
    
    ### Select the best model based on BIC
    best_full_index <- which.min(bic_full_list)
    fit_list_full[[c]][[i]] <- model_full_list[[best_full_index]]
    
    ### Run each algorithm N times and compute BIC for partial model
    for (n in 1:N) {
      fit_c_partial = tryCatch(
        fitMLCWMd(data = data, formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
                  C = c, 
                  cont_var =  c("x1", "x2"), 
                  cat_var = c("V1", "D1", "D2", "D3"), 
                  dich_dep_var = NULL, 
                  init = "kmeans", 
                  mstep_isingfit = TRUE, ###
                  fixed_intercept = TRUE,
                  max_it = 30, 
                  tol = 1e-6),
        error = function(e) {
          list(model = NA)
        }
      )
      
      # Calculate BIC if model is not NA
      if (is.null(fit_c_partial$model)) {
        bic_partial_list[n] <- fit_c_partial$bic
      } else {
        bic_partial_list[n] <- Inf  # Set to infinity if model fails
      }
      
      model_partial_list[[n]] <- fit_c_partial
    }
    
    ### Select the best model based on BIC
    best_partial_index <- which.min(bic_partial_list)
    fit_list_partial[[c]][[i]] <- model_partial_list[[best_partial_index]]
    
  }
  
  ### Save results
  saveRDS(fit_list_full,"fitted_models/fit_list_full")
  saveRDS(fit_list_partial,"fitted_models/fit_list_partial")
  
  ### Check
  print(i)
}


#---------------------Get comparison models results-----------------------------
fit_glm <- list()
fit_glmer <- list()
for(i in 1:R){
  
  set.seed(i*7)
  data = readRDS( paste0("sim_data/mydata", i, ".rds") )
  
  fit_glm[[i]] = tryCatch(
    glm(formula = y ~ x1 + x2 + V1 + D1 + D2 + D3, data = data, 
        family = binomial),
    error = function(e) {
      list(model = NA)
    }
  )
  
  fit_glmer[[i]] = tryCatch(
    glmer(data = data, formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
          family = binomial),
    error = function(e) {
      list(model = NA)
    }
  )
  print(i)
}
saveRDS(fit_glm,"fitted_models/fit_list_glm")
saveRDS(fit_glmer,"fitted_models/fit_list_glmer")
