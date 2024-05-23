# Model Prediction

model_pred <- function(fit, U, V = NULL, D = NULL, C, data, new_level = FALSE){
  v = length(V)
  weights_ = matrix(0, nrow = dim(data)[1], ncol = C)
  for(i in 1:C){
    weights_[, i] = exp( log(fit$params$w[i]) + dmvnorm(U, mean = 
                                                          fit$params$mu[,,i], sigma = fit$params$sigma[,,i],log = 1) +
                           + all_mydmultinom_log(V,fit$params$lambda,v,i)  +
                           log(apply(D, 1, function(row) IsingStateProb(row, 
                                                                                 fit$params$int[,,i], fit$params$thres[,,i], beta = 1, 
                                                                                 responses = c(0L, 1L)))) ) 
    
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


table_Pred <- function(y_pred, data_pred){
  N <- length(y_pred)
  border_string <- "border-bottom: 2px solid;border-top: 2px solid;border-left: 2px solid;border-right: 2px solid;"
  ff = as.matrix( cbind(data_pred, Prediction = round(y_pred,2)) )
  cn <- c("New_Obs", colnames(ff) )
  colnames(ff) <- NULL
  ff %>% kbl(align = c("c", "c","c", "c","c", "c","c", "c","c", "c","c","c", "c","c", "c","c")) %>%
    kable_styling(font_size = 20)  %>%
    row_spec(seq(1,length(y_pred)), extra_css = "border-bottom: 2px solid;") %>%
    column_spec(1:16, extra_css = "border-left: 2px solid;") %>%
    add_header_above(cn, extra_css = border_string, 
                     background = alpha("lightblue",0.2))
  
}

