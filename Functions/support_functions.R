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


# cat <--> one-hot

onehot_to_cat <- function(data, one_hot_cols_list, new_names) {
  if (length(one_hot_cols_list) != length(new_names)) {
    stop("The length of one_hot_cols_list must match the length of new_names.")
  }
  
  for (i in seq_along(one_hot_cols_list)) {
    one_hot_cols <- one_hot_cols_list[[i]]
    new_name <- new_names[i]
    
    if (!all(one_hot_cols %in% names(data))) {
      stop(paste("One or more specified one-hot encoded columns for", new_name, "do not exist in the data."))
    }
    
    data[[new_name]] <- colnames(data[one_hot_cols])[max.col(data[one_hot_cols] == 1, ties.method = "first")]
    data <- data[, !(names(data) %in% one_hot_cols)]
  }
  
  return(data)
}


cat_to_onehot <- function(data, categorical_vars, new_names_list) {
  if (length(categorical_vars) != length(new_names_list)) {
    stop("The length of categorical_vars must match the length of new_names_list.")
  }
  
  for (i in seq_along(categorical_vars)) {
    categorical_var <- categorical_vars[i]
    new_names <- new_names_list[[i]]
    
    if (!categorical_var %in% names(data)) {
      stop(paste("The specified categorical variable", categorical_var, "does not exist in the data."))
    }
    
    data[[categorical_var]] <- as.factor(data[[categorical_var]])
    existing_levels <- levels(data[[categorical_var]])
    if (length(new_names) != length(existing_levels)) {
      stop(paste("The length of new_names for", categorical_var, "must match the total number of unique levels."))
    }
    one_hot <- model.matrix(~ get(categorical_var) - 1, data = data)
    colnames(one_hot) <- new_names
    data <- cbind(data, one_hot)
    data <- data[, !(names(data) %in% categorical_var)]
  }
  
  return(data)
}
