# Model Prediction
predict_MLCWMd <- function(fit, data, C, new_level = FALSE){
  if(!is.null(fit$params$mu)){
     U <- data[ , dimnames(fit$params$mu)[[2]] ]
     cont_term <- list()
     for(i in 1:C){
         cont_term[[i]] <- dmvnorm(U, mean = 
                                     fit$params$mu[,,i], sigma = fit$params$sigma[,,i],log = 1)
       }
  }  else{
    cont_term <- rep(0,C)
  }
  
  if(!is.null(fit$params$lambda)){
    V <- list()
    v <- length(fit$params$lambda)
    for(i in 1:v){
      temp_var <- sub("\\..*", "", dimnames(fit$params$lambda[[i]])[[2]])[1]
      temp_data <- as.data.frame(data[, temp_var])
      colnames(temp_data) <- temp_var
      V[[i]] <- cat_to_onehot(temp_data, temp_var, list(dimnames(fit$params$lambda[[i]])[[2]]))
    }
    multinom_term <- list()
    for(i in 1:C){
        multinom_term[[i]] <- all_mydmultinom_log(V,fit$params$lambda,v,i) 
      } 
  } else{
    multinom_term <- rep(0,C)
  }
  
  if(!is.null(fit$params$thres)){
     D <- data[ , dimnames(fit$params$thres)[[2]]]
     dich_term <- list()
     for(i in 1:C){
         dich_term[[i]] <- log(apply(D, 1, function(row) IsingStateProb(row, 
                  fit$params$int[,,i], fit$params$thres[,,i], beta = 1, responses = c(0L, 1L))))
     }
   } else{
         dich_term <- rep(0,C)
   } 
  
  
  weights_ = matrix(0, nrow = dim(data)[1], ncol = C)
  for(i in 1:C){
       weights_[, i] = exp( log(fit$params$w[i]) + cont_term[[i]] +
                              + multinom_term[[i]] + dich_term[[i]] ) 
       
  }
  y_pred <- 0
  if(new_level == F){
    for(i in 1:C){
      y_pred <- y_pred + weights_[,i]*predict(fit$params$models[[i]], data, type = "response", allow.new.levels = F) 
      
    }
  }
  else{
    for(i in 1:C){
      y_pred <- y_pred +  weights_[,i]*predict(fit$params$models[[i]], data, type = "response", allow.new.levels = T) 
      
    }
  }
  
  y_pred <- y_pred / rowSums(weights_)
  
  return(y_pred)
}

table_Pred <- function(y_pred, data_pred, fit) {
  if(!is.null(fit$params$mu)){
    u <- length(names(fitting$params$mu[,,1]))
  } else{
    u <- 0
  }
  if(!is.null(fit$params$thres)){
    d <- length(names(fitting$params$thres[,,1]))
  } else{
    d <- 0
  }
  if(!is.null(fit$params$lambda)){
    v <- length(fitting$params$lambda)
  } else{
    v <- 0
  }
  N <- length(y_pred)
  border_string <- "font-family: cursive; border-bottom: 2px solid;border-top: 2px solid;border-left: 2px solid;border-right: 2px solid;"
  ff = as.matrix( cbind("New" = 1:N, Prediction = round(y_pred,2), data_pred ) )
  cn <- colnames(ff)
  colnames(ff) <- NULL
  ff %>% kbl(align = c("c", "c","c", "c","c", "c","c", "c","c", "c","c","c", "c","c", "c","c")) %>%
      kable_styling(font_size = 20)  %>%
      row_spec(seq(1,length(y_pred)), extra_css = "border-bottom: 2px solid;") %>%
      column_spec(2, background = "#42B3FF") %>%
      column_spec(1, background = alpha("grey90",0.2)) %>%
      column_spec(1:(d + u + v + 3 + 1), extra_css = "border-left: 2px solid;", width = "100px") %>%
      column_spec(d + u + v + 3 + 1, extra_css = "border-right: 2px solid;") %>%
      add_header_above(cn, extra_css = border_string, color = "black", 
                       background = alpha("white",0.2))
}


