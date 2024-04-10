Algo_full <- function(data,C,U,V = NULL,D = NULL, flag_init = TRUE, flag_mstep = TRUE){
  
  # Decide Initialization type
  if(flag_init == T){
  params = initialize_params_random(data,U,V,D,formula,C)
  }
  else{
    params = initialize_params_k_means(data,U,V,D,formula,C)
  }
  
  z <- E_step(Y,U,V,D,C,params)
  log_l <- c(0,log_like(Y,U,V,D,C,z,params))
  tol = 1e-6
  while(abs(log_l[length(log_l)] - log_l[length(log_l)-1]) >  tol){
    # E-step
    z <- E_step(Y,U,V,D,C,params)
    # M-step
    if(flag_mstep == T){
         params <- M_step_EI(formula,U,V,D,data,C,z,params)
    }
    else{
         params <- M_step_IF(formula,U,V,D,data,C,z,params)
    }
    #Saving log-likelihood
    new_log_l <- log_like(Y,U,V,D,C,z,params)
    log_l <- c(log_l,new_log_l)
  }
  return(list(params = params,log_l = log_l,z = z))
}