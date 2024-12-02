fitMLCWMd <- function(data, formula, C, cont_var = NULL, cat_var = NULL, dich_dep_var = NULL, 
                      init, cluster = NULL, mstep_isingfit = FALSE, 
                      fixed_intercept = TRUE,
                      max_it = 30, tol = 1e-6){
  
  # Process variables
  if(length(cont_var) != 0){
    missing_cont_var <- setdiff(cont_var, colnames(data))
    if(length(missing_cont_var) == 0){
      U <- data[ , cont_var]
    } else{
      return("At least one continuous variable provided is not in the data")
    }
  } else{
    U <- NULL
  }
  
  if(length(cat_var) != 0){
    missing_cat_var <- setdiff(cat_var, colnames(data))
    if(length(missing_cat_var) == 0 ){
      cat_names_list <- list()
      for(i in 1:length(cat_var)){
        #cat_names_list[[i]] <- paste0(cat_var[i], 1:length( levels( as.factor(data[, cat_var[i]]) ) ) )
        cat_names_list[[i]] <- paste0(cat_var[i], ".", levels( as.factor(data[, cat_var[i]])  ) )
      }
      cat_data <- as.data.frame(data[ , cat_var])
      colnames(cat_data) <- cat_var
      cat_data <- as.data.frame(unclass(cat_data),stringsAsFactors = TRUE)
      onehot_cat_data <- cat_to_onehot(cat_data, cat_var, cat_names_list)
      V <- list()
      for(i in 1:length(cat_var)){
        V[[i]] <- onehot_cat_data[, grep(paste0("^", cat_var[i]), names(onehot_cat_data))]
      }
    } else{
      return("At least one categorical variable provided is not in the data")
    }
  } else{
    V <- NULL
  } 
  
  if(length(dich_dep_var) != 0){
     missing_dich_dep_var <- setdiff(dich_dep_var, colnames(data))
     if(length(missing_dich_dep_var) == 0){
       D <- data[ , dich_dep_var]
     } else{
       return("At least one dichotomous dependent variable provided is not in the data")
     }
  } else{
    D <- NULL
  }
  
  Y <- data[ , as.character(formula[[2]])]
  data <- as.data.frame(unclass(data),stringsAsFactors = TRUE)
  #return(list(U,V,D))
  
  # Decide Initialization type
  if(init == "random"){
    params = initialize_params_random(data,U,V,D,formula,C)
  } else if (init == "kmeans"){
    params = initialize_params_k_means(data,U,V,D,formula,C)
  } else if (init == "manual"){
    params = initialize_params_manual(data,cluster,U,V,D,formula,C)
  } else{
    stop("Execution stopped. No existing initialization provided")
  }
  
  z <- E_step(Y,U,V,D,C,params)
  log_l <- c(0,log_like(Y,U,V,D,C,z,params))
  tol = tol
  max_it = max_it
  n_iter = 0
  while( ( abs(log_l[length(log_l)] - log_l[length(log_l)-1])  >  tol ) & ( n_iter < max_it ) ){
    # E-step
    z <- E_step(Y,U,V,D,C,params)
    # M-step
    if(mstep_isingfit == FALSE){
      params <- M_step_EI(formula,U,V,D,data,C,z,params)
    }
    else{
      params <- M_step_IF(formula,U,V,D,data,C,z,params)
    }
    #Saving log-likelihood
    new_log_l <- log_like(Y,U,V,D,C,z,params)
    log_l <- c(log_l,new_log_l)
    n_iter <- n_iter + 1
  }
  model_bic = Calculate_BIC(formula, data, U, V, D, log_l[length(log_l)], C, fixed_intercept = fixed_intercept)
  
  if(length(cont_var) != 0){
      dimnames(params$mu) <- list(NULL, cont_var, NULL)
      dimnames(params$sigma) <- list(cont_var, cont_var, NULL)
  } 
  
  if(length(cat_var) != 0){
     for(i in 1:length(cat_var)){
        dimnames(params$lambda[[i]]) <- list(NULL, cat_names_list[[i]], NULL)
     }
  }
  if(length(dich_dep_var) != 0){
     dimnames(params$thres) <- list(NULL, dich_dep_var, NULL)
     dimnames(params$int) <- list(dich_dep_var, dich_dep_var, NULL)
  }
  
  return(list(params = params,log_l = log_l, z = z, bic = model_bic))
}
