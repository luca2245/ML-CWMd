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

fit_list_full = vector("list",4)
fit_list_partial = vector("list",4)

R = 100
dir.create("fitted_models")

for (i in 1:R){
  
  set.seed(i*7)
  
  data = readRDS( paste0("sim_data/mydata", i, ".rds") )
  
  for (c in 1:4){
    
    ### Fit with full model
    
    fit_c_full = tryCatch(
      fitMLCWMd(data = data, formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
                C = c, 
                cont_var =  c("x1", "x2"), 
                cat_var = c("V1"), 
                dich_dep_var = c("D1", "D2", "D3"), 
                init = "kmeans", 
                mstep_isingfit = TRUE, 
                fixed_intercept = TRUE,
                max_it = 30, 
                tol = 1e-6)
      ,
      error = function(e) {
        list(model = NA)
      }
    )
    
    fit_list_full[[c]][[i]] = fit_c_full
    
    ### Fit without dichotomous dependent
    
    fit_c_partial = tryCatch(
      fitMLCWMd(data = data, formula = y ~ x1 + x2 + V1 + D1 + D2 + D3 + (1|level),
                C = c, 
                cont_var =  c("x1", "x2"), 
                cat_var = c("V1", "D1", "D2", "D3"), 
                dich_dep_var = NULL, 
                init = "kmeans", 
                mstep_isingfit = TRUE, 
                fixed_intercept = TRUE,
                max_it = 30, 
                tol = 1e-6)
      ,
      error = function(e) {
        list(model = NA)
      }
    )
    
    fit_list_partial[[c]][[i]] = fit_c_partial
    
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
  saveRDS(fit_glm,"fitted_models/fit_list_glm")
  saveRDS(fit_glmer,"fitted_models/fit_list_glmer")
  
  print(i)
}
