
# Function for ML-CWMd BIC calculation
Calculate_BIC <- function(U, V = NULL, D = NULL, log_l, C, fixed_intercept = FALSE){
  N = dim(U)[1]
  p = dim(U)[2]
  l = dim(D)[2]
  k_r = unlist(purrr::map(V, ~dim(.)[2]))
  
  k_weights = C - 1
  k_cont = p*(p+3)/2
  k_cat = sum(k_r)
  k_dich = l*(l+1)/2
  k_reg = length(attr(terms(formula), "term.labels")) + ifelse(fixed_intercept == TRUE, 1, 0)
  
  k = C*(k_cont + k_cat + k_dich + k_reg) + k_weights
  
  BIC <- -2*log_l + log(N) * k
  return(BIC)
  
}


  
  
  