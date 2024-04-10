# E-step Function

E_step <- function(Y,U,V = NULL,D = NULL,C,params){
  list2env(params, env = environment())
  v <- length(V)
  if( length(D) != 0){
    dich_density <- list()
    for(i in 1:C){
      dich_density[[i]] <- log(apply(D, 1, function(row) IsingStateProb(row, 
                                                                        int[,,i], thres[,,i], beta = 1, responses = c(0L, 1L))))  
    }
  }
  else{
    dich_density <- rep(0,C)
  }
  if( length(V) != 0){
    multinom_density <- list()
    for(i in 1:C){
      multinom_density[[i]] <- all_mydmultinom_log(V,lambda,v,i)  
    }
  }
  else{
    multinom_density <- rep(0,C)
  }
  
  z = matrix(0, nrow = dim(U)[1], ncol = C)
  for(i in 1:C){
    z[, i] = exp( log(w[i]) + cross_log(fitted_values[,,i],Y) +
                    dmvnorm(U, mean = mu[,,i], sigma = sigma[,,i],log = 1) +
                    + multinom_density[[i]] +
                    dich_density[[i]]  )
  }
  z <- z / rowSums(z)
  z = mclust::unmap(mclust::map(z))
  return (z)
}