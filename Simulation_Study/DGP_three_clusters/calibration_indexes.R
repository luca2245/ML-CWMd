
# Expected Calibration Error
getECE <- function(actual, predicted, n_bins=10){ #equal frequency bins
  
  predicted <- predicted
  labels <- actual
  idx <- order(predicted)
  pred_actual <- (cbind(predicted[idx], labels[idx]))
  
  N <- nrow(pred_actual)
  rest <- N%%n_bins
  S <- 0
  W <- c()
  B <- min(N,n_bins) #if less then n_bins elements in data set, then use that number of bins
  groups <- list()
  
  for (i in 1:B){ #i von 1 bis B
    if (i <= rest){ #put rest elements into each bin
      group_pred <- (pred_actual[(((i-1)*ceiling(N/n_bins)+1) : (i*ceiling(N/n_bins))),1])
      group_actual <- (pred_actual[(((i-1)*ceiling(N/n_bins)+1) : (i*ceiling(N/n_bins))),2])
    }
    else {
      group_pred <- (pred_actual[((rest+(i-1)*floor(N/n_bins)+1) : (rest+i*floor(N/n_bins))),1])#group size=N/B
      group_actual <- (pred_actual[((rest+(i-1)*floor(N/n_bins)+1) : (rest+i*floor(N/n_bins))),2])
    }
    
    n_ <- length(group_pred)
    expected <- mean(group_pred) #mean of predictions in bin b
    observed <- mean(group_actual) #true fraction of pos.instances = prevalence in bin b
    
    S[i] <- abs(observed-expected) #absolut difference of observed value-predicted value in bin
    W[i] <- n_/N #empirical frequence of all instances that fall into bin i, should be equal when using equal freq binning approach
    groups[[i]] <- group_pred
    
  }
  
  mean_prediction <- lapply(groups, mean)
  min_group <- lapply(groups, min)
  max_group <- lapply(groups, max)
  
  res <- t(S)%*%W
  return(as.numeric(res))
}

# Maximum Calibration Error
getMCE <- function(actual, predicted, n_bins=10){
  
  predicted <- predicted
  labels <- actual
  idx <- order(predicted)
  pred_actual <- (cbind(predicted[idx], actual[idx]))
  N <- nrow(pred_actual)
  rest <- N%%n_bins
  B <- min(N,n_bins)
  
  S <- 0
  W <- c()
  for (i in 1:B){ #i von 1 bis B
    if (i <= rest){ #put rest elements into each bin
      group_pred <- (pred_actual[(((i-1)*ceiling(N/n_bins)+1) : (i*ceiling(N/n_bins))),1])
      group_actual <- (pred_actual[(((i-1)*ceiling(N/n_bins)+1) : (i*ceiling(N/n_bins))),2])
    }
    else {
      group_pred <- (pred_actual[((rest+(i-1)*floor(N/n_bins)+1) : (rest+i*floor(N/n_bins))),1])#group size=N/B
      group_actual <- (pred_actual[((rest+(i-1)*floor(N/n_bins)+1) : (rest+i*floor(N/n_bins))),2])
    }
    
    n <- length(group_pred)
    expected <- mean(group_pred) #mean of predictions in bin b
    observed <- mean(group_actual) #true fraction of pos.instances = prevalence in bin b
    
    S[i] <- abs(observed-expected) #absolut difference of observed value-predicted value in bin
    W[i] <- n/N #empirical frequence of all instances that fall into bin i, should be pretty much the same among all bins
  }
  
  res <- max(S*W)
  return(res)
}

# Brier Score
get_Brier_score <- function(actual, predicted){
  n <- length(actual)
  n_1 <- length(actual==1)
  n_0 <- length(actual==0)
  sum <- 0
  sum_0 <- 0
  sum_1 <- 0
  
  for (i in seq(1,n,1)){
    diff <- abs((predicted[i]-actual[i]))^2
    sum <- sum+diff
    if(actual[i]==0){
      sum_0 <- sum_0+diff
    }
    else
      if(actual[i]==1){
        sum_1 <- sum_1+diff
      }
  }
  return(list(brier=sum/n, brier_1=sum_1/n_1, brier_0=sum_0/n_0))
  
}


