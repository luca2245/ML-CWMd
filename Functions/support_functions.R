# SUPPORT FUNCTIONS

cross_log <- function(pi,y){
  interm <- ( pi^y * (1 - pi)^(1-y) )
  return (log(interm))
}


#Receive a list of covariates and a list of vectors of probabilities and return the 
#probability of observing each data point according to all covariates.

#- j: cluster
#- n: number of categorical covariates
#- V: list of covariates

all_mydmultinom_log <- function(V,prob,v,j){
  result = rep(0,dim(V[[1]])[1])
  for(i in 1:v){
    result = result + apply(V[[i]], 1, function(row) 
      dmultinom(x = row, size = 1, prob = prob[[i]][,,j], log = TRUE))
  }
  return(result)
}

# Multiplications for parameters update in the M-step:

mymult3 <- function(w,arr1,U){
  result <- matrix(0,nrow = dim(U)[2],ncol = dim(U)[2])
  for(i in 1:dim(arr1)[1]){
    result <- result + w[i] * arr1[i,] %*% t(arr1[i,])
  }
  return(result)
}



mymult4 <- function(w,V){
  result <- colSums(w*V)
  return(as.vector(result))
}
