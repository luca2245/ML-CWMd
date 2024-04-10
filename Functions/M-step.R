# M-step Function with EstimateIsing
M_step_EI <- function(formula,U,V = NULL,D = NULL,data,C,z,params){
  
  list2env(params, env = environment())
  new_params = list()
  sum_z <- colSums(z)
  
  #Update Mixture weights w
  w = sum_z / dim(data)[1]
  new_params$w <- w
  
  #Update Parameters related to U
  mu = t(as.matrix(z)) %*% as.matrix(U)
  mu = array(t(mu / sum_z), dim=c(1, dim(U)[2], C))
  for (c in 1:C){
    j = as.matrix(U - matrix(mu[,,c],nrow=dim(U)[1],ncol=dim(U)[2],byrow=T))
    s = mymult3(z[,c],j,U)
    s = s / sum_z[c]
    sigma[,,c] = s
  }
  new_params$mu <- mu
  new_params$sigma <- sigma
  
  #Update Parameters related to V
  if( length(V) != 0){
    v <- length(V)
    for(c in 1:C){
      for(j in 1:v){
        lambda[[j]][,,c] <- mymult4(z[,c],V[[j]])/sum_z[c]
      }
      
    }
    new_params$lambda <- lambda
  }
  
  #Update parameters related to D
  if( length(D) != 0){
    for (c in 1:C){
      ising_fit <- EstimateIsing(D[which(z[,c] == 1),],beta = 1,responses = c(0L, 1L), method = "pl")
      int[,,c] <- ising_fit$graph
      thres[,,c] <- ising_fit$thresholds
    }
    new_params$int <- int
    new_params$thres <- thres
  }
  
  #Update Parameters related to Y
  for (c in 1:C){
    result = do.call(glmer, 
                     list(formula=formula,data=data,family=binomial,weights=z[,c]))
    fitted_values[,,c] <- fitted(result)
    models[[c]] <- result
  }
  new_params$fitted_values <- fitted_values
  new_params$models <- models
  
  return(new_params)
}


# M-step Function with IsingFit
M_step_IF <- function(formula,U,V = NULL,D = NULL,data,C,z,params){
  
  list2env(params, env = environment())
  new_params = list()
  sum_z <- colSums(z)
  
  #Update Mixture weights w
  w = sum_z / dim(data)[1]
  new_params$w <- w
  
  #Update Parameters related to U
  mu = t(as.matrix(z)) %*% as.matrix(U)
  mu = array(t(mu / sum_z), dim=c(1, dim(U)[2], C))
  for (c in 1:C){
    j = as.matrix(U - matrix(mu[,,c],nrow=dim(U)[1],ncol=dim(U)[2],byrow=T))
    s = mymult3(z[,c],j,U)
    s = s / sum_z[c]
    sigma[,,c] = s
  }
  new_params$mu <- mu
  new_params$sigma <- sigma
  
  #Update Parameters related to V
  if( length(V) != 0){
    v <- length(V)
    for(c in 1:C){
      for(j in 1:v){
        lambda[[j]][,,c] <- mymult4(z[,c],V[[j]])/sum_z[c]
      }
      
    }
    new_params$lambda <- lambda
  }
  
  #Update parameters related to D
  if( length(D) != 0){
    for (c in 1:C){
      ising_fit <- IsingFit(D[which(z[,c] == 1),], family='binomial', plot=FALSE, progressbar = FALSE)
      int[,,c] <- ising_fit$weiadj
      thres[,,c] <- ising_fit$thresholds
    }
    new_params$int <- int
    new_params$thres <- thres
  }
  
  
  #Update Parameters related to Y
  for (c in 1:C){
    result = do.call(glmer, 
                     list(formula=formula,data=data,family=binomial,weights=z[,c]))
    fitted_values[,,c] <- fitted(result)
    models[[c]] <- result
  }
  new_params$fitted_values <- fitted_values
  new_params$models <- models
  
  return(new_params)
}