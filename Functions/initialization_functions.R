# K-MEANS Initialization
initialize_params_k_means <- function(data, U = NULL, V = NULL, D = NULL, formula, C){
  
  params = list()
  #k-means
  kmeans_result <- kmeans(U, centers = C)
  #Initialize w[c] for any c as the proportion of data points assigned to 
  #cluster c by k-means
  
  U <- data.frame(U)
  w <- vector(mode = "numeric", length = C)
  for(c in 1:C){
    w[c] <- sum(kmeans_result$cluster == c)/dim(U)[1]
  }
  params$w = w
  
  #Initialize mu and sigma as the centers and the covariance matrix returned by k-means
  if(length(U) != 0){
     mu <- array(0, dim=c(1, dim(U)[2], C))
     for(c in 1:C){
       mu[,,c] <- kmeans_result$centers[c,]
     }
     sigma <- array(0, dim=c(dim(U)[2], dim(U)[2], C))
     for(c in 1:C){
       if(dim(U)[2] > 1){
         sigma[,,c] <- cov(U[which(kmeans_result$cluster == c),])
       } else{
         sigma[,,c] <- var(U[which(kmeans_result$cluster == c),])
       }
      
     }
     params$mu <- mu
     params$sigma <- sigma
  }
  #Initialize each vector lambda_{vc} (c: c-th latent group, v: v-th covariate)
  #with information of k-means
  if(length(V) != 0){
    
    v <- length(V)
    lambda <- vector(mode = "list", length = v)
    for(i in 1:v){
      lambda[[i]] <- array(0, dim = c(1,dim(V[[i]])[2],C))
      for (c in 1:C) {
        vec <- V[[i]][which(kmeans_result$cluster == c),]
        lambda[[i]][,,c] <- colSums(vec)/dim(vec)[1]  
      }
    }
    params$lambda <- lambda
  }
  
  #Initialize parameters related to dependent dichotomic variables with the Ising Estimate in the groups
  #identified by initial_z
  initial_z <- mclust::unmap(kmeans_result$cluster)
  if(length(D) != 0){
    
    d = dim(D)[2]
    thres <- array(0, dim=c(1, d, C))
    int <- array(0, dim=c(d, d, C))
    for (c in 1:C){
      ising_fit <- IsingFit(D[which(initial_z[,c] == 1),], family='binomial', plot=FALSE, progressbar = FALSE)
      int[,,c] <- ising_fit$weiadj
      thres[,,c] <- ising_fit$thresholds
    }
    params$int <- int
    params$thres <- thres
  }
  
  #Initialize the parameters and fitted values of the multilevel logistic model 
  #in each cluster exploiting the information of k-means, in particular un-mapping the 
  #clusters' vector returned by k-means
  fitted_values <- array(0,dim = c(dim(data)[1],1,C))
  models <- vector(mode = "list", length = C)
  xi <- array(list(), dim = c(1, C))
  for (c in 1:C){
    result = do.call(glmer, 
                     list(formula=formula,data=data,family=binomial,weights=initial_z[,c]))
    fitted_values[,,c] <- fitted(result)
    models[[c]] <- result
  }
  params$fitted_values <- fitted_values
  params$models <- models
  
  return(params)
}

#RANDOM Initialization with Check
initialize_params_random <- function(data, U = NULL, V = NULL, D = NULL, formula, C){

  params = list()
  
  N <- dim(data)[1]
  random_z <- sample(C, N, replace = TRUE)
  check_vec <- vector(mode = "numeric", length = C)
  for(c in 1:C){
    check_vec[c] <- sum(random_z == c)/N
  }
  while(sum(check_vec > 1/(C+1) ) != C){
    random_z <- sample(C, N, replace = TRUE)
    for(c in 1:C){
       check_vec[c] <- sum(random_z == c)/N
  }
   }


  w <- vector(mode = "numeric", length = C)
  for(c in 1:C){
    w[c] <- sum(random_z == c)/N
  }
  params$w = w

  #Initialize mu and sigma as the centers and the covariance matrix returned by k-means
  U <- data.frame(U)
  if(length(U) != 0){
     mu <- array(0, dim=c(1, dim(U)[2], C))
     for(c in 1:C){
       mu[,,c] <- colMeans(U[which(random_z == c),])
     }
     sigma <- array(0, dim=c(dim(U)[2], dim(U)[2], C))
     for(c in 1:C){
       sigma[,,c] <- cov(U[which(random_z == c),])
     }
     params$mu <- mu
     params$sigma <- sigma
  }
  #Initialize each vector lambda_{vc} (c: c-th latent group, v: v-th covariate)
  #with information of k-means
  if(length(V) != 0){

    v <- length(V)
    lambda <- vector(mode = "list", length = v)
    for(i in 1:v){
      lambda[[i]] <- array(0, dim = c(1,dim(V[[i]])[2],C))
      for (c in 1:C) {
        vec <- V[[i]][which(random_z == c),]
        lambda[[i]][,,c] <- colSums(vec)/dim(vec)[1]
      }
    }
    params$lambda <- lambda
  }

  #Initialize parameters related to dependent dichotomic variables with the Ising Estimate in the groups
  #identified by initial_z
  if(length(D) != 0){
    d <- N <- dim(D)[2]
    int1 <- matrix(sample(0:1,N^2,TRUE,prob = c(0.5, 0.5)),N,N) * rnorm(N^2)
    int1 <- int1 + t(int1)
    diag(int1) <- 0
    thres1 <- rnorm(d,0,0.5)
    thres <- array(0, dim=c(1, d, C))
    int <- array(0, dim=c(d, d, C))
    for (c in 1:C){
      int[,,c] <- int1
      thres[,,c] <- thres1
    }
    params$int = int
    params$thres = thres
  }

  #Initialize the parameters and fitted values of the multilevel logistic model
  #in each cluster exploiting the information of k-means, in particular un-mapping the
  #clusters' vector returned by k-means
  models <- vector(mode = "list", length = C)
  lmm1 = glmer(formula,family = binomial,
               data = data,
               glmerControl(boundary.tol = 1e-2,check.conv.grad = .makeCC("warning",
                                                                          tol = 7e-2, relTol = NULL)) )
  fitted_values <- array(c(rep(fitted(lmm1),C)),dim = c(dim(data)[1],1,C))
  for (c in 1:C){
    models[[c]] <- lmm1
  }
  params$fitted_values = fitted_values
  params$models = models


  return(params)
}


#Manual Initialization
initialize_params_manual <- function(data, cluster, U = NULL, V = NULL, D = NULL, formula, C){
  params = list()
  
  N <- dim(data)[1]
  #Initialize w[c] for any c as the proportion of data points assigned to
  #cluster c by k-means
  w <- vector(mode = "numeric", length = C)
  for(c in 1:C){
    w[c] <- sum(cluster == c)/N
  }
  params$w = w

  #Initialize mu and sigma as the centers and the covariance matrix returned by k-means
  if(length(U) != 0){
    mu <- array(0, dim=c(1, dim(U)[2], C))
    for(c in 1:C){
      mu[,,c] <- mean(U[which(cluster == c),])
    }
    sigma <- array(0, dim=c(dim(U)[2], dim(U)[2], C))
    for(c in 1:C){
      sigma[,,c] <- cov(U[which(cluster == c),])
    }
    params$mu <- mu
    params$sigma <- sigma
  }
  #Initialize each vector lambda_{vc} (c: c-th latent group, v: v-th covariate)
  #with information of k-means
  if(length(V) != 0){

    v <- length(V)
    lambda <- vector(mode = "list", length = v)
    for(i in 1:v){
      lambda[[i]] <- array(0, dim = c(1,dim(V[[i]])[2],C))
      for (c in 1:C) {
        vec <- V[[i]][which(cluster == c),]
        lambda[[i]][,,c] <- colSums(vec)/dim(vec)[1]
      }
    }
    params$lambda <- lambda
  }

  #Initialize parameters related to dependent dichotomic variables with the Ising Estimate in the groups
  #identified by initial_z
  initial_z <- mclust::unmap(cluster)
  if(length(D) != 0){

    d = dim(D)[2]
    thres <- array(0, dim=c(1, d, C))
    int <- array(0, dim=c(d, d, C))
    for (c in 1:C){
      ising_fit <- IsingFit(D[which(initial_z[,c] == 1),], family='binomial', plot=FALSE, progressbar = FALSE)
      int[,,c] <- ising_fit$weiadj
      thres[,,c] <- ising_fit$thresholds
    }
    params$int <- int
    params$thres <- thres
  }

  #Initialize the parameters and fitted values of the multilevel logistic model
  #in each cluster exploiting the information of k-means, in particular un-mapping the
  #clusters' vector returned by k-means
  fitted_values <- array(0,dim = c(dim(data)[1],1,C))
  models <- vector(mode = "list", length = C)
  xi <- array(list(), dim = c(1, C))
  for (c in 1:C){
    result = do.call(glmer,
                     list(formula=formula,data=data,family=binomial,weights=initial_z[,c]))
    fitted_values[,,c] <- fitted(result)
    models[[c]] <- result
  }
  params$fitted_values <- fitted_values
  params$models <- models

  return(params)
}

