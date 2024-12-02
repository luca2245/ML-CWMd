
# Function for ML-CWMd BIC calculation
Calculate_BIC <- function(formula, data, U = NULL, V = NULL, D = NULL, log_l, C, fixed_intercept = FALSE){
  N = dim(data)[1]
  k_weights = C - 1
  
  if(!is.null(U)){
    U <- data.frame(U)
    p = dim(U)[2]
    k_cont = p*(p+3)/2
  } else{
    k_cont <- 0
  }
  
  if(!is.null(V)){
    k_r = unlist(purrr::map(V, ~dim(.)[2]))
    k_cat = sum(k_r)
  } else{
    k_cat <- 0
  }
  
  if(!is.null(D)){
    l = dim(D)[2]
    k_dich = l*(l+1)/2
  } else{
    k_dich <- 0
  }
  
  k_reg = length(attr(terms(formula), "term.labels")) + ifelse(fixed_intercept == TRUE, 1, 0)
  
  k = C*(k_cont + k_cat + k_dich + k_reg) + k_weights
  
  BIC <- -2*log_l + log(N) * k
  return(BIC)
  
}


  
  
  