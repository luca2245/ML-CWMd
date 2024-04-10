# Log-Likelihood Function
log_like <- function(Y,U,V = NULL,D = NULL,C,z,params){
  
  list2env(params, env = environment())
  v <- length(V)
  d = dim(D)[2]
  
  if( length(D) != 0){
    dich_Lik <- list()
    for(i in 1:C){
      dich_Lik[[i]] <- log(apply(D, 1, function(row) IsingStateProb(row, 
                                                                    int[,,i], thres[,,i], beta = 1, responses = c(0L, 1L))))  
    }
  }
  else{
    dich_Lik <- rep(0,C)
  }
  if( length(V) != 0){
    multinom_Lik <- list()
    for(i in 1:C){
      multinom_Lik[[i]] <- all_mydmultinom_log(V,lambda,v,i)  
    }
  }
  else{
    multinom_Lik <- rep(0,C)
  }
  
  tot = 0
  for (c in 1:C){
    tot = tot +  z[,c] * ( log(w[c])  + 
                             dmvnorm(U, mean = mu[,,c], sigma = sigma[,,c],log = TRUE) +
                             multinom_Lik[[c]] + dich_Lik[[c]] )
  }
  tot_lik_glmer = 0
  for (c in 1:C){
    tot_lik_glmer <- tot_lik_glmer + logLik(models[[c]])[1]
  }
  
  
  
  return ( sum(tot) + tot_lik_glmer)
}