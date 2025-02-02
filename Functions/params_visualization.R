library(kableExtra)
library(gridExtra)

pedix <- c("\u2081","\u2082","\u2083","\u2084","\u2085","\u2086","\u2087","\u2088","\u2089")
border_string <- "border-bottom: 2px solid;border-top: 2px solid;border-left: 2px solid;border-right: 2px solid;"
my_colors <- c("rgba(0, 0, 255, 0.6)", "rgba(0, 255, 0, 0.6)",
               "rgba(255, 20, 147, 0.6)", "rgba(255, 165, 0, 0.6)", "rgba(255, 0, 0, 0.6)",
               "rgba(160, 32, 240, 0.6)")

# Process data
process_data <- function(data,
                         cont_var = c("x1", "x2"), 
                         cat_var = c("V1"), 
                         dich_dep_var = c("D1", "D2", "D3")){
  # Process variables
  if(length(cont_var) != 0){
    missing_cont_var <- setdiff(cont_var, colnames(data))
    if(length(missing_cont_var) == 0){
      U <- data[ , cont_var]
    } else{
      return("At least one continuous variable provided is not in the data")
    }
  } else{
    U <- NULL
  }
  
  if(length(cat_var) != 0){
    missing_cat_var <- setdiff(cat_var, colnames(data))
    if(length(missing_cat_var) == 0 ){
      cat_names_list <- list()
      for(i in 1:length(cat_var)){
        #cat_names_list[[i]] <- paste0(cat_var[i], 1:length( levels( as.factor(data[, cat_var[i]]) ) ) )
        cat_names_list[[i]] <- paste0(cat_var[i], ".", levels( as.factor(data[, cat_var[i]])  ) )
      }
      cat_data <- as.data.frame(data[ , cat_var])
      colnames(cat_data) <- cat_var
      cat_data <- as.data.frame(unclass(cat_data),stringsAsFactors = TRUE)
      onehot_cat_data <- cat_to_onehot(cat_data, cat_var, cat_names_list)
      V <- list()
      for(i in 1:length(cat_var)){
        V[[i]] <- onehot_cat_data[, grep(paste0("^", cat_var[i]), names(onehot_cat_data))]
      }
    } else{
      return("At least one categorical variable provided is not in the data")
    }
  } else{
    V <- NULL
  } 
  
  if(length(dich_dep_var) != 0){
    missing_dich_dep_var <- setdiff(dich_dep_var, colnames(data))
    if(length(missing_dich_dep_var) == 0){
      D <- data[ , dich_dep_var]
    } else{
      return("At least one dichotomous dependent variable provided is not in the data")
    }
  } else{
    D <- NULL
  }
  
  return(list(U = U, V = V, D = D))
}




# MU params
visual_mu_params <- function(fitting){
   mu <- data.frame(matrix(ncol = 0,nrow = dim(U)[2]))
   mu$Variable <- names(fitting$params$mu[,,1])
   for(c in 1:C){
     coln <- paste0("mu","_",c)
     mu[[coln]] = signif(fitting$params$mu[,,c], digits = 4)
   }
   colnames(mu) <- c("Variable", paste0("\u03BC",pedix[1:C]) ) 
   mu1 <- as.matrix(mu)
   colnames(mu1) <- NULL
   kable_mu <- mu1 %>% kbl(align = c("c","c","c","c")) %>%
     kable_styling(font_size = 20) %>% 
     column_spec(1, color = "black",
                 background = alpha("grey",0.2), bold = TRUE)  %>% 
     column_spec(2, color = "black",
                 background = my_colors[1]) %>%
     column_spec(3, color = "black",
                 background = my_colors[2]) 
   
   if(C == 3 | C == 4 | C == 5 | C == 6){
     kable_mu <- kable_mu %>% column_spec(4, color = "black",
                                          background = my_colors[3])
   }
   if(C == 4 | C == 5 | C == 6){
     kable_mu <- kable_mu %>% column_spec(5, color = "black",
                                          background = my_colors[4])
   }
   if(C == 5 | C == 6){
     kable_mu <- kable_mu %>% column_spec(6, color = "black",
                                          background = my_colors[5])
   }
   if(C == 6){
     kable_mu <- kable_mu %>% column_spec(7, color = "black",
                                          background = my_colors[6])
   }
   
   kable_mu <- kable_mu %>% column_spec(1:(C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
     row_spec(seq(1,dim(U)[2]), extra_css = "border-bottom: 2px solid;") %>%
     column_spec(1, extra_css = "border-left: 2px solid;") %>%
     add_header_above(colnames(mu), extra_css = border_string) 
   return(kable_mu)
}

# SIGMA params
visual_sigma_params <- function(fitting){
  sigma <- data.frame(matrix(ncol = 0,nrow = dim(U)[2]))
  var_names <- names(fitting$params$mu[,,1])
  sigma$Variable <- var_names
  for(c in 1:C){
    sigma <- cbind(sigma,signif(fitting$params$sigma[,,c],digits = 4),
                   names(fitting$params$mu[,,c]))
  }
  sigma <- sigma[,-ncol(sigma)]
  sig_vec <- paste0("\u03A3",pedix[1:C])
  vec <- NULL
  for(i in 1:C){
    vec <- c(vec,"Variables", rep(sig_vec[i],dim(U)[2]) )
  }
  vec2 <- NULL
  for(i in 1:C){
    vec2 <- c(vec2," ", rep(sig_vec[i],dim(U)[2]) )
  }
  colnames(sigma) <- vec
  tt <- data.frame(matrix(  rep(c("Variables",
                                  var_names),C) ,ncol = dim(U)[2]*C + C,
                            nrow = 1)  )
  colnames(tt) <- vec
  sigma2 <- rbind(tt,sigma)
  sigma1 <- as.matrix(sigma2)
  colnames(sigma2) <- NULL
  rownames(sigma2) <- NULL
  step_size <- dim(U)[2] - 1
  kable_sigma <- sigma2 %>% kbl(align = c("c","c","c","c")) %>%
    kable_styling(font_size = 20) %>% column_spec(seq(1,dim(U)[2]*C + C,dim(U)[2]+1), color = "black",
                                                  background = alpha("grey",0.2), bold = TRUE)  %>%
    column_spec(2:(2+step_size), color = "black",
                background = my_colors[1]) %>%
    column_spec((4+step_size):(4+2*step_size), color = "black",
                background = my_colors[2])
  
  if(C == 3 | C == 4 | C == 5 | C == 6){
    kable_sigma <- kable_sigma %>% column_spec((6+2*step_size):(6+3*step_size),
                                               color = "black",
                                               background = my_colors[3])
  }
  if(C == 4 | C == 5 | C == 6){
    kable_sigma <- kable_sigma %>% column_spec((8+3*step_size):(8+4*step_size), color = "black",
                                               background = my_colors[4])
  }
  if(C == 5 | C == 6){
    kable_sigma <- kable_sigma %>% column_spec((10+4*step_size):(10+5*step_size), color = "black",
                                               background = my_colors[5])
  }
  if(C == 6){
    kable_sigma <- kable_sigma %>% column_spec((12+5*step_size):(12+6*step_size), color = "black",
                                               background = my_colors[6])
  }
  
  kable_sigma <- kable_sigma %>% row_spec(1, color = "black",
                                          background = alpha("grey",0.2), bold = TRUE) %>%
    column_spec(seq(1,dim(U)[2]*C + C,dim(U)[2]+1), border_right = "2px solid black", border_left = "2px solid black") %>%
    column_spec(dim(U)[2]*C + C, border_right = "2px solid black") %>%
    row_spec(c(1,dim(U)[2]+1), extra_css = "border-bottom: 2px solid;") %>%
    add_header_above(vec2)
  return(kable_sigma)
}


# LAMBDA params
visual_lambda_params <- function(fitting){
  V1 <- V[[1]]
  fir <- data.frame(matrix(ncol = 0,nrow = dim(V1)[2]))
  fir$Category <- sub(".*\\.", "", colnames(V1))
  fir1 <- data.frame(cbind(lambda_11 =  signif(fitting$params$lambda[[1]][,,1],digits = 4),
                           lambda_12 = signif(fitting$params$lambda[[1]][,,2],digits = 4)) )
  if(C == 3 | C == 4 | C == 5 | C == 6){
    fir1 <- cbind(fir1, lambda_13 = signif(fitting$params$lambda[[1]][,,3],digits = 4))
  }
  if(C == 4 | C == 5 | C == 6){
    fir1 <- cbind(fir1, lambda_14 = signif(fitting$params$lambda[[1]][,,4],digits = 4))
  }
  if(C == 5 | C == 6){
    fir1 <- cbind(fir1, lambda_15 = signif(fitting$params$lambda[[1]][,,5],digits = 4))
  }
  if(C == 6){
    fir1 <- cbind(fir1, lambda_16 = signif(fitting$params$lambda[[1]][,,6],digits = 4))
  }
  
  fir <- cbind(fir,fir1)
  colnames(fir) <- c("Category", paste0("\u03BB",pedix[1], pedix[1:C]) )
  fir1 <- as.matrix(fir)
  colnames(fir1) <- NULL
  rownames(fir1) <- NULL
  kable_fir1 <- fir1 %>% kbl(align = c("c","c","c","c","c")) %>%
    kable_styling(font_size = 20) %>% 
    column_spec(1, color = "black",
                background = alpha("grey",0.2)) %>%
    column_spec(2, color = "black",
                background = my_colors[1]) %>%
    column_spec(3, color = "black",
                background = my_colors[2])
  
  if(C == 3 | C == 4  | C == 5 | C == 6){
    kable_fir1 <- kable_fir1  %>% column_spec(4, color = "black",
                                              background = my_colors[3])
  }
  if(C == 4| C == 5 | C == 6){
    kable_fir1 <- kable_fir1  %>% column_spec(5, color = "black",
                                              background = my_colors[4])
  }
  if(C == 5 | C == 6){
    kable_fir1 <- kable_fir1  %>% column_spec(6, color = "black",
                                              background = my_colors[5])
  }
  if(C == 6){
    kable_fir1 <- kable_fir1  %>% column_spec(7, color = "black",
                                              background = my_colors[6])
  }
  
  kable_fir1 <- kable_fir1 %>%
    column_spec(seq(1,C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
    row_spec(seq(1,dim(V1)[2]), extra_css = "border-bottom: 2px solid;") %>%
    add_header_above(colnames(fir), extra_css = border_string) 
  
  return(kable_fir1)
}


# NU params
visual_nu_params <- function(fitting){
  nuu <- data.frame(matrix(ncol = 0,nrow = dim(D)[2]  ))
  nuu$Variable <- paste0("D",1:dim(D)[2])
  for(c in 1:C){
    coln <- paste0("nu","\u208",c)
    nuu[[coln]] = signif(fitting$params$thres[,,c], digits = 3)
  }
  colnames(nuu) <- c("Variable", paste0("\u03BD",pedix[1:C]) )
  nuu1 <- as.matrix(nuu)
  colnames(nuu1) <- NULL
  kable_nuu <- nuu1 %>% kbl(align = c("c","c","c")) %>%
    kable_styling(font_size = 20) %>% column_spec(1, color = "black",
                                                  background = alpha("grey",0.2)) %>%
    column_spec(2, color = "black",
                background = my_colors[1]) %>%
    column_spec(3, color = "black",
                background = my_colors[2]) 
  
  if(C == 3 | C == 4 | C == 5 | C == 6){
    kable_nuu <- kable_nuu %>% column_spec(4, color = "black",
                                           background = my_colors[3])
  }
  if(C == 4 | C == 5 | C == 6){
    kable_nuu <- kable_nuu %>% column_spec(5, color = "black",
                                           background = my_colors[4])
  }
  if(C == 5 | C == 6){
    kable_nuu <- kable_nuu %>% column_spec(6, color = "black",
                                           background = my_colors[5])
  }
  if(C == 6){
    kable_nuu <- kable_nuu %>% column_spec(7, color = "black",
                                           background = my_colors[6])
  }
  
  kable_nuu <- kable_nuu %>%
    column_spec(1:(C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
    row_spec(seq(1,dim(D)[2]), extra_css = "border-bottom: 2px solid;") %>%
    add_header_above(colnames(nuu), extra_css = border_string) 
  
  return(kable_nuu)
}


# THR params
visual_thr_params <- function(fitting) {
  lab <- paste0("D", 1:dim(fitting$params$int)[1])  
  plots <- list()  
  my_colors <- c("blue", "green", "deeppink", "orange", "red", "purple")
  for (i in 1:C) {
    qgraph(fitting$params$int[,,i], fade = FALSE, 
           edge.labels = round(fitting$params$int[,,i], 3), 
           edge.label.cex = 1.5, labels = lab, edge.label.color = "black",
           color = alpha(my_colors[i], 0.5),  # Assign dynamic colors
           edge.width = 0.7, negCol = "red", posCol = "forestgreen", 
           title = paste0("Cluster ", i, " Fitted network"), 
           title.cex = 2.5, edge.label.margin = 0.001, 
           label.font = 2, edge.label.font = 2)
    
    plots[[i]] <- recordPlot()
  }
  
  return(plots)  
}


# FIXED EFFECTS
visual_fixef_params <- function(fitting) {
  tables <- list()
  colors <- c("rgba(0, 0, 255, 0.6)", "rgba(0, 255, 0, 0.6)",
              "rgba(255, 20, 147, 0.6)", "rgba(255, 165, 0, 0.6)", "rgba(255, 0, 0, 0.6)",
              "rgba(160, 32, 240, 0.6)")
  
  for (i in 1:C) {
    # Extract coefficients and format them
    fixef_data <- as.data.frame(signif(summary(fitting$params$models[[i]])$coefficients, digits = 3))
    fixef_matrix <- as.matrix(fixef_data)
    colnames(fixef_matrix) <- NULL  # Remove column names for display
    
    # Create the table using kableExtra
    table <- fixef_matrix %>%
      kbl(align = c("c", "c", "c", "c", "c", "c")) %>%
      kable_styling(font_size = 15) %>% 
      row_spec(which(fixef_data[, 4] < 0.05), color = "black", background = colors[i]) %>%
      row_spec(which(fixef_data[, 4] >= 0.05), color = "black", background = alpha("grey", 0.3)) %>%
      column_spec(1, color = "black", bold = TRUE) %>%
      column_spec(1:5, border_right = "2px solid black", border_left = "2px solid black") %>%
      row_spec(seq(1, nrow(fixef_data)), extra_css = "border-bottom: 2px solid;") %>%
      add_header_above(c("Variables", colnames(fixef_data)), extra_css = border_string)
    
    tables[[i]] <- table  # Store the table in the list
  }
  
  return(tables)  # Return the list of tables
}

# RANDOM EFFECTS
visual_ranef_params <- function(fitting) {
  
  
  plots <- list()
  
  my_colors <- c("blue", "green", "deeppink", "orange", "red", "purple")
  for (i in seq_along(fitting$params$models)) {
    
    re_sim <- REsim(fitting$params$models[[i]])
    re_sim$CI_lower <- re_sim$mean - (1.96 * re_sim$sd)
    re_sim$CI_upper <- re_sim$mean + (1.96 * re_sim$sd)
    re_sim$contains_zero <- with(re_sim, CI_lower <= 0 & CI_upper >= 0)
    re_sim$groupFctr <- " "
    re_sim <- re_sim %>% arrange(median)
    
    # Set colors based on confidence interval containing zero
    x_colors <- ifelse(re_sim$contains_zero, "grey30", my_colors[i])
    
    # Generate plot
    plot <- plotREsim(re_sim, labs = TRUE) +
      geom_hline(yintercept = 0, color = I("black"), size = 1.2) +
      geom_errorbar(aes(color = contains_zero), width = 0, size = 0.8) +
      geom_point(aes(color = contains_zero), size = 3) +
      scale_color_manual(values = c(my_colors[i], "gray50"), guide = "none") +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 16, face = "bold", vjust = 3, angle = 0, color = x_colors),
        plot.title = element_text(color = my_colors[i], size = 20, vjust = -2),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 16, face = "bold", hjust = 1.5)
      )
    
    # Store plot in the list
    plots[[i]] <- plot
  }
  
  return(plots)
}







