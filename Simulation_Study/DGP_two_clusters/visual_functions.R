library(reshape2)
library(ggtext)
library(clue) 
library(progress)
library(CalibrationCurves)
source("calibration_indexes.R")
#---------------------------------Clusters plot---------------------------------
cluster_plot <- function(data, fit, add_to_title = "",
                         my_colors = c("#6666FF", "#CC0066", "#009966", "#FF3300", "#660099")){
  my_colors_vec <- c("1" = my_colors[1], "2" = my_colors[2],  "3" = my_colors[3],  "4" = my_colors[4], 
                     "5" = my_colors[5])
  p <- ggplot(data) +
    geom_point(size = 4, aes(x = x1,
                             y = x2, 
                             fill = as.factor(mclust::map(fit$z))), shape = 21) +
    scale_fill_manual(values = alpha(my_colors_vec, 0.6)) +
    labs(x = names(fit$params$mu[,,1])[1],
         y = names(fit$params$mu[,,1])[2],
         fill = "Cluster", title = paste0("Clusters Plot", add_to_title)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_bw() +
    theme(plot.title = element_text(size = 24, face = "bold.italic", hjust = 0.5),
          axis.title = element_text(size = 24, face = "bold.italic"),
          legend.position = "right",
          legend.title = element_text(face = "bold.italic", size = 20),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 11) ) +
    guides(color = "none")
  
  return(p)
}

cluster_plot_byclust <- function(data, clusters, add_to_title = "",
                         my_colors = c("#6666FF", "#CC0066", "#009966", "#FF3300", "#660099")){
  my_colors_vec <- c("1" = my_colors[1], "2" = my_colors[2],  "3" = my_colors[3],  "4" = my_colors[4], 
                     "5" = my_colors[5])
  p <- ggplot(data) +
    geom_point(size = 4, aes(x = x1,
                             y = x2, 
                             fill = as.factor(clusters)), shape = 21) +
    scale_fill_manual(values = alpha(my_colors_vec, 0.6)) +
    labs(x = "x1",
         y = "x2",
         fill = "Cluster", title = paste0("Clusters Plot", add_to_title)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_bw() +
    theme(plot.title = element_text(size = 24, face = "bold.italic", hjust = 0.5),
          axis.title = element_text(size = 24, face = "bold.italic"),
          legend.position = "right",
          legend.title = element_text(face = "bold.italic", size = 20),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 11) ) +
    guides(color = "none")
  
  return(p)
}

#------------------------------ML-CWMd comparison-------------------------------
MLCWMd_comparison <- function(full, partial, glm_res = NULL, glmer_res = NULL, x_breaks, title, y_title,
                              my_colors = c("forestgreen", "darkorange", "red", "blue"), 
                              col_x_axis = FALSE, show.outliers = T){
  x_labb <- c("ML-CWMd", "ML-CWMd \nwith no D")
  j <- 2
  data <- data.frame(cbind("ML-CWMd" = full, "ML-CWMd \nwith no D" =  partial))
  if(!is.null(glmer_res)){
    x_labb <- c(x_labb, "GLMM")
    data <- cbind(data, "GLMM" = glmer_res)
    j <- j + 1
  }
  if(!is.null(glm_res)){
    x_labb <- c(x_labb, "GLM")
    data <- cbind(data, "GLM" = glm_res)
    j <- j + 1
  }
  my_colors <- my_colors[1:j]
  melted_data <- melt(data, id.vars = NULL)
  
  my_colors <- my_colors
  p1 <- ggplot(melted_data , aes(x = variable, y = value, fill = variable, color = variable)) +
    geom_boxplot(outliers = show.outliers) +
    scale_fill_manual(values = alpha(my_colors, 0.3)) +  # Use adjusted colors
    scale_color_manual(values = alpha(my_colors, 1)) +
    labs(title = title,
         x = "Model",
         y = y_title,
         fill = "Model") +
    scale_x_discrete(labels = x_labb) +
    scale_y_continuous(breaks = x_breaks) +
    theme_minimal() +  
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold.italic"),
      axis.title = element_text(size = 12, face = "bold.italic"),
      legend.position = "none",
      legend.text = element_text(size = 12),
      legend.title = element_text(face = "bold.italic", size = 14),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 11) 
    ) +
    guides(color = "none") 
  
  if(col_x_axis == T){
    p1 <- p1 + theme(axis.text.x = ggtext::element_markdown(size = 12, face = "bold.italic", colour = my_colors) )
    
  } else{
    p1 <- p1 +  theme(axis.text.x = element_text(size = 12, face = "bold.italic"))  
  }
  return(p1)
}

calculate_misclassification <- function(fit_list, sim_data_list, C_real, n_sim) {
  miss <- vector(mode = "numeric", length = n_sim)
  
  for (k in 1:n_sim) {
    # Check if model exists
    if (is.null(fit_list[[C_real]][[k]]$model)) {
      fitted_clusters <- mclust::map(fit_list[[C_real]][[k]]$z)
      true_clusters <- sim_data_list[[k]]$latent
      
      tab_obs <- table(fitted_clusters, true_clusters)
      perm <- solve_LSAP(tab_obs, maximum = TRUE) 
      
      matched_confusion <- tab_obs[cbind(1:nrow(tab_obs), perm)]
      misclassified <- sum(tab_obs) - sum(matched_confusion) 
      miss[k] <- misclassified / sum(tab_obs)  
    } else {
      miss[k] <- NA  
    }
  }
  
  return(miss)
}


BIC_comparison_plot <- function(bic_list, n_sim, lwd,
                    my_colors = c("#6666FF", "#CC0066", "#009966", "#FF3300"),
                    lty = "solid", v_just_ann = -0.2, h_just_ann = -1.8){
  my_colors_lighter <- c(alpha(my_colors, 0.5), alpha(my_colors, 0.2) )
  my_colors_ <- rep(my_colors, 2)
  max_length <- max(sapply(bic_list, length))
  bic_list_padded <- lapply(bic_list, function(x) {
    length(x) <- max_length
    return(x)
  })
  data <- as.data.frame(bic_list_padded)
  colnames(data) <- as.character(1:length(bic_list))
  
  melted_data <- melt(data, id.vars = NULL)
  cl_ass <- rep(1:length(bic_list), each = n_sim)
  data_for_plot <- cbind(melted_data, cl_ass)
  
  p <- ggplot(data_for_plot, aes(x = variable, y = value, fill = as.factor(cl_ass), 
                                 color = as.factor(cl_ass) ) ) +
    geom_boxplot(width = 0.8) 
  
  breaks_col <- NULL
  labels_col <- NULL
  values_col <- NULL
  values_fill <- NULL
  for(i in 1:length(bic_list)){
    breaks_col <- c(breaks_col, as.character(i) )
    labels_col <- c(labels_col, names(bic_list)[i] )
    values_col <- c(values_col, setNames(my_colors_[i], as.character(i) ))
    values_fill <- c(values_fill, setNames(my_colors_lighter[i], as.character(i) ))
  }
  p <- p + scale_color_manual(name = 'Model', breaks = breaks_col, 
                              labels = labels_col, values = values_col )
  p <- p + scale_fill_manual(name = 'Model', breaks = breaks_col, 
                             labels = labels_col, values = values_fill )
  p <- p +
    labs(title = "BIC Comparison",
         x = "Parameters",
         y = "BIC", fill = "Model", color= "Model") +
    scale_x_discrete(labels = names(bic_list)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24, face = "bold.italic", hjust = 0.5),
      axis.title.y = element_text(size = 20, face = "bold.italic"),
      axis.title.x = element_blank(),
      legend.position = "none",
      legend.text = element_text(size = 18),
      legend.key.width = unit(3, "lines"),
      legend.title = element_text(face = "bold.italic", size = 18),
      axis.text.x = element_text(size = 18, face = "bold.italic"),
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 11) 
    ) +
    guides(color = "none")
  
  return(p)
}

calculate_ARI <- function(fit_list, sim_data_list, C_real, n_sim) {
  ari <- vector(mode = "numeric", length = n_sim)
  for(k in 1:n_sim){
    if(is.null(fit_list[[C_real]][[k]]$model)){
      ari[k] <- adjustedRandIndex(mclust::map(fit_list[[C_real]][[k]]$z), 
                                  sim_data_list[[k]]$latent)
    } else{
      ari[k] <- NA
    }
  } 
  return(ari)
}

get_MLCWMd_pred_accuracy <- function(fit_list, sim_data_list, test_data, C_real, n_sim) {
  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed",
    total = n_sim,
    clear = FALSE,
    width = 80
  )
  acc_train <- vector(mode = "numeric", length = n_sim)
  acc_test <- vector(mode = "numeric", length = n_sim)
  for(k in 1:n_sim){
    if(is.null(fit_list[[C_real]][[k]]$model)){
      train_pred <- predict_MLCWMd(fit_list[[C_real]][[k]], sim_data_list[[k]], C_real, new_level = FALSE)
      test_pred <- predict_MLCWMd(fit_list[[C_real]][[k]], test_data, C_real, new_level = FALSE)
      roc_curve <- roc(sim_data_list[[k]][,"y"], train_pred, levels = c(0, 1),
                       direction = "<")
      optimal.threshold <- coords(roc_curve, "best", best.method = "closest.topleft")
      
      predicted_classes_train <- ifelse(train_pred > as.numeric(optimal.threshold[1]), 1, 0)
      acc_train[k] <- mean(predicted_classes_train == sim_data_list[[k]][,"y"])
      
      predicted_classes_test <- ifelse(test_pred > as.numeric(optimal.threshold[1]), 1, 0)
      acc_test[k] <- mean(predicted_classes_test == test_data[,"y"])
    } else{
      acc_train[k] <- NA
      acc_test[k] <- NA
    }
    
    pb$tick()
  } 
  return(list(train_acc = acc_train, test_acc = acc_test))
}


get_glm_pred_acc <- function(fit_list, sim_data_list, test_data, n_sim) {
  acc_train <- vector(mode = "numeric", length = n_sim)
  acc_test <- vector(mode = "numeric", length = n_sim)
  for(k in 1:n_sim){
    predicted_probs_train <- predict(fit_list[[k]], type = "response")
    roc_curve <- roc(sim_data_list[[k]][,"y"], predicted_probs_train, levels = c(0, 1),
                     direction = "<")
    optimal.threshold <- coords(roc_curve, "best", best.method = "closest.topleft")
    predicted_classes_train <- ifelse(predicted_probs_train > as.numeric(optimal.threshold[1]), 1, 0)
    acc_train[k] <- mean(predicted_classes_train == sim_data_list[[k]][,"y"])
    
    predicted_probs_test <- predict(fit_list[[k]], test_data, type = "response")
    predicted_classes_test <- ifelse(predicted_probs_test > as.numeric(optimal.threshold[1]), 1, 0)
    acc_test[k] <- mean(predicted_classes_test == test_data[,"y"])
  }
  return(list(train_acc = acc_train, test_acc = acc_test))
}

get_glmer_pred_acc <- function(fit_list, sim_data_list, test_data, n_sim) {
  acc_train <- vector(mode = "numeric", length = n_sim)
  acc_test <- vector(mode = "numeric", length = n_sim)
  for(k in 1:n_sim){
    predicted_probs_train <- predict(fit_list[[k]], type = "response")
    roc_curve <- roc(sim_data_list[[k]][,"y"], predicted_probs_train, levels = c(0, 1),
                     direction = "<")
    optimal.threshold <- coords(roc_curve, "best", best.method = "closest.topleft")
    predicted_classes_train <- ifelse(predicted_probs_train > as.numeric(optimal.threshold[1]), 1, 0)
    acc_train[k] <- mean(predicted_classes_train == sim_data_list[[k]][,"y"])
    
    predicted_probs_test <- predict(fit_list[[k]], test_data, type = "response")
    predicted_classes_test <- ifelse(predicted_probs_test > as.numeric(optimal.threshold[1]), 1, 0)
    acc_test[k] <- mean(predicted_classes_test == test_data[,"y"])
  }
  return(list(train_acc = acc_train, test_acc = acc_test))
}


# Get Calibration
get_glm_calib <- function(fit_list, test_data, n_sim, flag = T) {
  ECE_idx <- vector(mode = "numeric", length = n_sim)
  MCE_idx <- vector(mode = "numeric", length = n_sim)
  ECI_idx <- vector(mode = "numeric", length = n_sim)
  Brier_idx <- vector(mode = "numeric", length = n_sim)
  if(flag == T){
    for(k in 1:n_sim){
      test_pred <- predict(fit_list[[k]], test_data, type = "response")
      calibcurv <- val.prob.ci.2(test_pred, test_data$y)
      ECE_idx[k] <- calibcurv$stats["Eavg"]
      MCE_idx[k] <- calibcurv$stats["Emax"]
      ECI_idx[k] <- calibcurv$stats["ECI"]
      Brier_idx[k] <- calibcurv$stats["Brier"]
    }
    return(list(ECE = ECE_idx, MCE = MCE_idx, ECI = ECI_idx, Brier = Brier_idx))
  } else{
    for(k in 1:n_sim){
      test_pred <- predict(fit_list[[k]], test_data, type = "response")
      ECE_idx[k] <- getECE(test_data$y, test_pred)
      MCE_idx[k] <- getMCE(test_data$y, test_pred)
      Brier_idx[k] <- get_Brier_score(test_data$y, test_pred)$brier
    }
    return(list(ECE = ECE_idx, MCE = MCE_idx, Brier = Brier_idx))
  }
}

get_glmer_calib <- function(fit_list, test_data, n_sim, flag = T) {
  ECE_idx <- vector(mode = "numeric", length = n_sim)
  MCE_idx <- vector(mode = "numeric", length = n_sim)
  ECI_idx <- vector(mode = "numeric", length = n_sim)
  Brier_idx <- vector(mode = "numeric", length = n_sim)
  if(flag == T){
    for(k in 1:n_sim){
      test_pred <- predict(fit_list[[k]], test_data, type = "response")
      calibcurv <- val.prob.ci.2(test_pred, test_data$y)
      ECE_idx[k] <- calibcurv$stats["Eavg"]
      MCE_idx[k] <- calibcurv$stats["Emax"]
      ECI_idx[k] <- calibcurv$stats["ECI"]
      Brier_idx[k] <- calibcurv$stats["Brier"]
    }
    return(list(ECE = ECE_idx, MCE = MCE_idx, ECI = ECI_idx, Brier = Brier_idx))
  } else{
    for(k in 1:n_sim){
      test_pred <- predict(fit_list[[k]], test_data, type = "response")
      ECE_idx[k] <- getECE(test_data$y, test_pred)
      MCE_idx[k] <- getMCE(test_data$y, test_pred)
      Brier_idx[k] <- get_Brier_score(test_data$y, test_pred)$brier
    }
    return(list(ECE = ECE_idx, MCE = MCE_idx, Brier = Brier_idx))
  }
}

get_MLCWMd_calib <- function(fit_list, test_data, C_real, n_sim, flag = T) {
  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed",
    total = n_sim,
    clear = FALSE,
    width = 80
  )
  ECE_idx <- vector(mode = "numeric", length = n_sim)
  MCE_idx <- vector(mode = "numeric", length = n_sim)
  ECI_idx <- vector(mode = "numeric", length = n_sim)
  Brier_idx <- vector(mode = "numeric", length = n_sim)
  if(flag == T){
    for(k in 1:n_sim){
      if(is.null(fit_list[[C_real]][[k]]$model)){
        test_pred <- predict_MLCWMd(fit_list[[C_real]][[k]], test_data, C_real, new_level = FALSE)
        calibcurv <- val.prob.ci.2(test_pred, test_data$y)
        ECE_idx[k] <- calibcurv$stats["Eavg"]
        MCE_idx[k] <- calibcurv$stats["Emax"]
        ECI_idx[k] <- calibcurv$stats["ECI"]
        Brier_idx[k] <- calibcurv$stats["Brier"]
      } else{
        ECE_idx[k] <- NA
        MCE_idx[k] <- NA
        ECI_idx[k] <- NA
        Brier_idx[k] <- NA
      }
      pb$tick()
    }
    return(list(ECE = ECE_idx, MCE = MCE_idx, ECI = ECI_idx, Brier = Brier_idx))
  } else{
    for(k in 1:n_sim){
      if(is.null(fit_list[[C_real]][[k]]$model)){
        test_pred <- predict_MLCWMd(fit_list[[C_real]][[k]], test_data, C_real, new_level = FALSE)
        ECE_idx[k] <- getECE(test_data$y, test_pred)
        MCE_idx[k] <- getMCE(test_data$y, test_pred)
        Brier_idx[k] <- get_Brier_score(test_data$y, test_pred)$brier
      } else{
        ECE_idx[k] <- NA
        MCE_idx[k] <- NA
        ECI_idx[k] <- NA
        Brier_idx[k] <- NA
      }
      pb$tick()
    }
    return(list(ECE = ECE_idx, MCE = MCE_idx, Brier = Brier_idx))
  }
}

#-----------Extract and plot real vs fitted parameters--------------------------
load_real_params <- function() {
  
  # Set DPG parameters
  n_obs <- 1200
  n_groups <- 8 # number of groups
  n_per_group <- 150 # number of observations per group
  group_ids <- rep(1:n_groups, each = n_per_group)
  w <- c(0.45,0.55) # proportion of observations per cluster
  
  # Random Intercept
  group_intercepts_1 <- c(0.69, -1.11, 0.17, -1.27, 1.14, -0.21, -2.83, 0.18)
  group_intercepts_2 <- c(0.57, 3.22, 0.86, 0.16, 0.16, -1.58, 1.11, 0.84)
  
  # MU and SIGMA Parameters
  cov.1 <- matrix(c(1,0.5,0.5,1),nrow = 2)
  mean.1 <- c(1.5, 0)
  
  cov.2 <- matrix(c(1,0.5,0.5,1),nrow = 2)
  mean.2 <- c(2.5, -0.8)
  
  # LAMBDA Parameters
  cat_probs1 <- c(0.25, 0.55, 0.2)
  cat_probs2 <- c(0.15, 0.65, 0.2)
  
  #INTERACTION AND THRESHOLD Parameters
  N <- 3
  
  int1 <- matrix(c(0, 1.8, 2.1, 1.8, 0, -2.4, 2.1, -2.4, 0), nrow = 3, ncol = 3)
  thres1 <- c(0, 0, 0)
  
  int2 <- matrix(c(0, -1.5, -3.0, -1.5,  0, 2.3, -3.0,  2.3, 0), nrow = 3, ncol = 3)
  thres2 <- c(0, 0, 0)
  
  # Fixed effects
  beta1 <- c(-0.6, -1.2, -0.8, 1.1, 0.9, 1.0, -0.7) 
  beta2 <- c(0.2, 0.4, 0.7, -0.6, -0.3, -1.3, 0.4)
  
  # Create a list of all variables to export
  params <- list(n_obs = n_obs, 
                 n_groups = n_groups, 
                 n_per_group = n_per_group, 
                 group_ids = group_ids, 
                 w = w,
                 group_intercepts_1 = group_intercepts_1, 
                 group_intercepts_2 = group_intercepts_2, 
                 cov.1 = cov.1, mean.1 = mean.1, 
                 cov.2 = cov.2, mean.2 = mean.2, 
                 cat_probs1 = cat_probs1, 
                 cat_probs2 = cat_probs2, 
                 N = N, 
                 int1 = int1, 
                 thres1 = thres1, 
                 int2 = int2, 
                 thres2 = thres2, 
                 beta1 = beta1, 
                 beta2 = beta2)
  
  # Assign all variables in the list to the global environment
  list2env(params, envir = .GlobalEnv)
}

match_clusters <- function(fitted_clust, true_clust){
  tab_obs <- table(true_clust, fitted_clust)
  perm <- solve_LSAP(tab_obs, maximum = TRUE) 
  perm <- cbind(1:nrow(tab_obs), perm)
  named_perm <- setNames(perm[, 2], perm[, 1])
  return(named_perm)
}

# Extracting continuous covariates params from simulated models
extract_mu_from_sim_data <- function(fit_list, sim_data_list, C){
  mu_dim <- dim(fit_list[[C]][[1]]$params$mu)[2]
  MU_list <- list()
  for(c in 1:C){
    MU_list[[c]] <- data.frame(matrix(0,nrow = length(sim_data_list), ncol = mu_dim))
    for(i in 1:length(sim_data_list)){
      if(is.null(fit_list[[C]][[i]]$model)){
        perm_vec <- match_clusters(mclust::map(fit_list[[C]][[i]]$z), sim_data_list[[i]]$latent)
        MU_list[[c]][i,] <- fit_list[[C]][[i]]$params$mu[, , perm_vec[as.character(c)]]  
      } else{
        MU_list[[c]][i,] <- NA
      }
    }
  }
  mu_data <- do.call(cbind, MU_list)
  col_names <- vector(mode = "character", length = mu_dim * C)
  i = 1
  for (c in 1:C) {
    for (j in 1:mu_dim) {
      col_names[i] <- paste0("mu_", c, j)
      i = i + 1
    }
  }
  colnames(mu_data) <- col_names
  return(mu_data)
}

extract_beta_from_sim_data <- function(fit_list, sim_data_list, C){
  beta_dim <- length(fit_list[[C]][[1]]$params$models[[1]]@beta) - 1
  BETA_list <- list()
  for(c in 1:C){
    BETA_list[[c]] <- data.frame(matrix(0,nrow = length(sim_data_list), ncol = beta_dim))
    for(i in 1:length(sim_data_list)){
      if(is.null(fit_list[[C]][[i]]$model)){
        perm_vec <- match_clusters(mclust::map(fit_list[[C]][[i]]$z), sim_data_list[[i]]$latent)
        BETA_list[[c]][i,] <- fit_list[[C]][[i]]$params$models[[perm_vec[as.character(c)]]]@beta[2:(beta_dim+1)] 
      } else{
        BETA_list[[c]][i,] <- NA
      }
    }
  }
  beta_data <- do.call(cbind, BETA_list)
  col_names <- vector(mode = "character", length = beta_dim * C)
  i = 1
  for (c in 1:C) {
    for (j in 1:beta_dim) {
      col_names[i] <- paste0("beta_", c, j)
      i = i + 1
    }
  }
  colnames(beta_data) <- col_names
  return(beta_data)
}


extract_sigma_from_sim_data <- function(fit_list, sim_data_list, C){
  sig_dim <- (dim(fit_list[[C]][[1]]$params$sigma)[1]*dim(fit_list[[C]][[1]]$params$sigma)[2]) - 1
  COV_list <- list()
  for(c in 1:C){
    COV_list[[c]] <- data.frame(matrix(0,nrow = length(sim_data_list), ncol = sig_dim))
    for(i in 1:length(sim_data_list)){
      if(is.null(fit_list[[C]][[i]]$model)){
        perm_vec <- match_clusters(mclust::map(fit_list[[C]][[i]]$z), sim_data_list[[i]]$latent)
        COV_list[[c]][i,] <- fit_list[[C]][[i]]$params$sigma[, , perm_vec[as.character(c)]][-3] 
      }
      else{
        COV_list[[c]][i,] <- NA
        
      }
    }
  }
  sigma_data <- do.call(cbind, COV_list)
  col_names <- vector(mode = "character", length = sig_dim * C)
  i = 1
  u <- dim(fit_list[[C]][[1]]$params$mu)[2]
  m_vec <- upper_diag_strings(u)
  for (c in 1:C) {
    for (j in 1:sig_dim) {
      col_names[i] <- paste0("Sigma", c, "_", m_vec[j])
      i = i + 1
    }
  }
  colnames(sigma_data) <- col_names
  return(sigma_data)
}

# Extract categorical covariates parameters
extract_lambda_from_sim_data <- function(fit_list, sim_data_list, C, v){
  lam_dim <- dim(fit_list[[C]][[1]]$params$lambda[[v]])[2]
  LAM_list <- list()
  for(c in 1:C){
    LAM_list[[c]] <- data.frame(matrix(0,nrow = length(sim_data_list), ncol = lam_dim))
    for(i in 1:length(sim_data_list)){
      if(is.null(fit_list[[C]][[i]]$model)){
        perm_vec <- match_clusters(mclust::map(fit_list[[C]][[i]]$z), sim_data_list[[i]]$latent)
        LAM_list[[c]][i,] <- fit_list[[C]][[i]]$params$lambda[[v]][, , perm_vec[as.character(c)]]
      }
      else{
        LAM_list[[c]][i,] <- NA
      }
    }
  }
  lam_data <- do.call(cbind, LAM_list)
  col_names <- vector(mode = "character", length = C * lam_dim)
  k = 1
  for (c in 1:C) {
    for (j in 1:lam_dim) {
      col_names[k] <- paste0("lambda_", c, v, j)
      k = k + 1
    }
  }
  colnames(lam_data) <- col_names
  return(lam_data)
}

# Extract dichotomous dependent covariates parameters
extract_inter_from_sim_data <- function(fit_list, sim_data_list, C, real_params_vec){
  N <- dim(fit_list[[C]][[1]]$params$int)[1]
  int_dim <- N*(N-1)/2
  INT_list <- list()
  for(c in 1:C){
    INT_list[[c]] <- data.frame(matrix(0,nrow = length(sim_data_list), ncol = int_dim))
    for(i in 1:length(sim_data_list)){
      if(is.null(fit_list[[C]][[i]]$model)){
        perm_vec <- match_clusters(mclust::map(fit_list[[C]][[i]]$z), sim_data_list[[i]]$latent)
        INT_list[[c]][i,] <- fit_list[[C]][[i]]$params$int[,, perm_vec[as.character(c)]][upper.tri(fit_list[[C]][[i]]$params$int[,,perm_vec[as.character(c)]])]
      }
      else{
        INT_list[[c]][i,] <- NA
      }
    }
  }
  int_data <- do.call(cbind, INT_list)
  col_names <- vector(mode = "character", length = C * int_dim)
  i = 1
  m_vec <- upper_diag_strings(int_dim)
  for (c in 1:C) {
    for (j in 1:int_dim) {
      col_names[i] <- paste0("gamma", c, "_", m_vec[j])
      i = i + 1
    }
  }
  colnames(int_data) <- col_names
  return(int_data)
}

extract_thres_from_sim_data <- function(fit_list, sim_data_list, C, real_params_vec){
  thres_dim <- dim(fit_list[[C]][[1]]$params$thres)[2]
  THR_list <- list()
  for(c in 1:C){
    THR_list[[c]] <- data.frame(matrix(0,nrow = length(sim_data_list), ncol = thres_dim))
    for(i in 1:length(sim_data_list)){
      if(is.null(fit_list[[C]][[i]]$model)){
        perm_vec <- match_clusters(mclust::map(fit_list[[C]][[i]]$z), sim_data_list[[i]]$latent)
        THR_list[[c]][i,] <- fit_list[[C]][[i]]$params$thres[, , perm_vec[as.character(c)]]
        
      } else{
        THR_list[[c]][i,] <- NA
      }
    }
  }
  thres_data <- do.call(cbind, THR_list)
  col_names <- vector(mode = "character", length = thres_dim * C)
  i = 1
  for (c in 1:C) {
    for (j in 1:thres_dim) {
      col_names[i] <- paste0("thres_", c, j)
      i = i + 1
    }
  }
  colnames(thres_data) <- col_names
  return(thres_data)
}

beta_plot <- function(data, comp_data, n_sim, lwd, C, vars,
                      my_colors = c("#6666FF", "#CC0066", "#009966", "#FF3300", "#660099"),
                      lty = "solid", v_just_ann = -0.2, h_just_ann = -1.8,
                      rm.outliers = TRUE){
  labels <- get_beta_labels(variables = vars)
  u <- length(vars)
  x_labels <- labels$x_labels     # Returns bold(mu[c*x[var]])
  x_labb_hat <- labels$x_labb_hat    # Returns hat(mu)[bold(c*x[var])]
  x_labb_real <- labels$x_labb_real   # Returns bold(mu[c*x[var]]^real)
  
  subdatasets <- split_dataset_by_columns(data, u, colnames(data)[1:u])
  melted_data <- NULL
  for (i in seq_along(subdatasets)) {
    melted <- melt(subdatasets[[i]], id.vars = NULL)
    melted_data <- rbind(melted_data, melted)
  }
  colnames(comp_data) <- colnames(data)[1:u] 
  melted_comp <- melt(comp_data, id.vars = NULL)
  
  x <- seq(0.55, by = 1, length.out = 1*u)
  xend <- seq(1.45, by = 1, length.out = 1*u)
  y <- yend <- do.call(c, lapply(1:C_real, function(i) get(paste0("beta", i))  ))
  cl_ass <- rep(1:C, each = n_sim * u)
  comp_mod_names <- c("GLM", "GLMM")
  comp_mod_ass <- rep(rep(comp_mod_names, each = n_sim), u)
  cl_ass <- c(cl_ass, comp_mod_ass)
  
  melted_data <- rbind(melted_data, melted_comp)
  data_for_plot <- cbind(melted_data, cl_ass)
  
  p <- ggplot(data_for_plot, aes(x = variable, y = value, fill = as.factor(cl_ass), color = as.factor(cl_ass) ) ) +
    geom_boxplot(aes(x = variable), outliers = rm.outliers, 
                 position = position_dodge(width = 0.9)) 
  
  breaks_col <- NULL
  labels_col <- NULL
  values_col <- NULL
  for(i in 1:C){
    breaks_col <- c(breaks_col, as.character(i) )
    labels_col <- c(labels_col, as.character(i) )
    values_col <- c(values_col, setNames(my_colors[i], as.character(i) ))
  }
  for(i in seq_along(comp_mod_names)){
    breaks_col <- c(breaks_col, comp_mod_names[i] )
    labels_col <- c(labels_col, comp_mod_names[i] )
    values_col <- c(values_col, setNames(my_colors[C + i], comp_mod_names[i] ))
  }
  ass_ <- rep(c(1:C), each = u)
  p <- p + purrr::map(1:u, function(i) {
    purrr::map(0:(C-1), function(j) {
      geom_segment(aes(x = x[i], y = y[i + u * j], xend = xend[i], yend = yend[i + u * j], 
                       color = as.character(ass_[i + u * j])), linetype = lty, linewidth = lwd)
    })
  })
  
  p <- p + scale_color_manual(name = 'Cluster', breaks = breaks_col, 
                              labels = labels_col, values = alpha(values_col, 0.7))
  p <- p + scale_fill_manual(name = 'Cluster', breaks = breaks_col, 
                             labels = labels_col, values = alpha(values_col, 0.3))
  p <- p +
    labs(subtitle = "Beta Parameters",
         title = "Regression Fixed Effects",
         x = "Parameters",
         y = "Estimated Values", fill = "Cluster", color= "Cluster") +
    scale_x_discrete(labels = x_labb_hat) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      plot.subtitle = element_text(size = 16, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold.italic"),
      legend.position = "right",
      legend.text = element_text(size = 18),
      legend.key.width = unit(3, "lines"),
      legend.title = element_text(face = "bold.italic", size = 18),
      axis.text.x = element_text(size = 24, face = "bold.italic"),
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 11) 
    ) +
    guides(color = "none")
  
  return(p)
}

mu_plot <- function(data, n_sim, lwd, C, vars,
                    my_colors = c("#6666FF", "#CC0066", "#009966", "#FF3300", "#660099"),
                    lty = "solid", v_just_ann = -0.2, h_just_ann = -1.8){
  labels <- get_mu_labels(C = C, variables = vars)
  u <- length(vars)
  x_labels <- labels$x_labels      # Returns bold(mu[c*x[var]])
  x_labb_hat <- labels$x_labb_hat    # Returns hat(mu)[bold(c*x[var])]
  x_labb_real <- labels$x_labb_real   # Returns bold(mu[c*x[var]]^real)
  
  x <- seq(0.6, by = 1, length.out = C*u)
  xend <- seq(1.4, by = 1, length.out = C*u)
  y <- yend <- do.call(c, lapply(1:C, function(i) get(paste0("mean.", i))))
  
  melted_data <- melt(data, id.vars = NULL)
  cl_ass <- rep(1:C, each = n_sim * u)
  data_for_plot <- cbind(melted_data, cl_ass)
  
  p <- ggplot(data_for_plot, aes(x = variable, y = value, fill = as.factor(cl_ass), color = as.factor(cl_ass) ) ) +
    geom_boxplot(width = 0.5) 
  
  breaks_col <- NULL
  labels_col <- NULL
  values_col <- NULL
  for(i in 1:C){
    breaks_col <- c(breaks_col, as.character(i) )
    labels_col <- c(labels_col, as.character(i) )
    values_col <- c(values_col, setNames(my_colors[i], as.character(i) ))
  }
  ass_ <- rep(c(1:C), each = u)
  p <- p + purrr::map(1:(C*u), function(i) {
    geom_segment(aes(x = x[i], y = y[i], xend = xend[i], yend = yend[i], 
                     color = as.character(ass_[i])), linetype = lty, linewidth = lwd)
     } )
  
  p <- p + scale_color_manual(name = 'Cluster', breaks = breaks_col, 
                              labels = labels_col, values = alpha(values_col, 0.7))
  p <- p + scale_fill_manual(name = 'Cluster', breaks = breaks_col, 
                              labels = labels_col, values = alpha(values_col, 0.3))
  p <- p +
    labs(subtitle = "Mean parameters",
         title = "Continuous covariates",
         x = "Parameters",
         y = "Estimated Values", fill = "Cluster", color= "Cluster") +
    scale_x_discrete(labels = x_labb_hat) +
    annotate("text", x = seq(1, by = 1, length.out = C*u), y = y, 
             label = x_labb_real, color = c(rep(my_colors[1:C], each = u)), 
             vjust = v_just_ann, hjust = h_just_ann, size = 5, parse = TRUE) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      plot.subtitle = element_text(size = 16, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold.italic"),
      legend.position = "right",
      legend.text = element_text(size = 18),
      legend.key.width = unit(3, "lines"),
      legend.title = element_text(face = "bold.italic", size = 18),
      axis.text.x = element_text(size = 22, face = "bold.italic", 
                                 color = rep(my_colors[1:C], each = 2) ),
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 11) 
    ) +
    guides(color = "none")
  
  return(p)
}


sigma_plot <- function(data, n_sim, lwd, C, vars,
              my_colors = c("#6666FF", "#CC0066", "#009966", "#FF3300", "#660099"),
              lty = "solid", v_just_ann = 0, h_just_ann = -1.0){
  u <- length(vars)
  labels <- get_sigma_labels(C = C, dims = u)
  
  x_labels <- labels$x_labels      # Returns bold(Sigma[cij])
  x_labb_hat <- labels$x_labb_hat    # Returns hat(Sigma)[bold(cij)]
  x_labb_real <- labels$x_labb_real   # Returns bold(Sigma[cij]^real)
  
  x <- seq(0.6, by = 1, length.out = C*(u*(u+1)/2) )
  xend <- seq(1.4, by = 1, length.out = C*(u*(u+1)/2) )
  y <- yend <- do.call(c, lapply(1:C_real, function(i) get(paste0("cov.", i))[-3] ))
  
  melted_data <- melt(data, id.vars = NULL)
  cl_ass <- rep(1:C, each = n_sim * u*(u+1)/2)
  data_for_plot <- cbind(melted_data, cl_ass)
  
  p <- ggplot(data_for_plot, aes(x = variable, y = value, fill = as.factor(cl_ass), color = as.factor(cl_ass) ) ) +
    geom_boxplot(width = 0.5) 
  
  breaks_col <- NULL
  labels_col <- NULL
  values_col <- NULL
  for(i in 1:C){
    breaks_col <- c(breaks_col, as.character(i) )
    labels_col <- c(labels_col, as.character(i) )
    values_col <- c(values_col, setNames(my_colors[i], as.character(i) ))
  }
  ass_ <- rep(c(1:C), each = u*(u+1)/2)
  p <- p + purrr::map(1:(C*u*(u+1)/2), function(i) {
    geom_segment(aes(x = x[i], y = y[i], xend = xend[i], yend = yend[i], 
                     color = as.character(ass_[i])), linetype = lty, linewidth = lwd)
  } )
  
  p <- p + scale_color_manual(name = 'Cluster', breaks = breaks_col, 
                              labels = labels_col, values = alpha(values_col, 0.7))
  p <- p + scale_fill_manual(name = 'Cluster', breaks = breaks_col, 
                             labels = labels_col, values = alpha(values_col, 0.3))
  
  p <- p + labs(subtitle = "Covariance matrix parameters",
                title = "Continuous covariates",
         x = "Parameters",
         y = "Estimated Values", fill = "Cluster", color= "Cluster") +
    scale_x_discrete(labels = x_labb_hat) +
    annotate("text", x = seq(1, by = 1, length.out = C*(u*(u+1)/2) ), y = y, 
             label = x_labb_real, color = c(rep(my_colors[1:C], each = u*(u+1)/2)), 
             vjust = v_just_ann, hjust = h_just_ann, size = 5, parse = TRUE) +
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      plot.subtitle = element_text(size = 16, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold.italic"),
      legend.position = "right",
      legend.text = element_text(size = 18),
      legend.key.width = unit(3, "lines"),
      legend.title = element_text(face = "bold.italic", size = 18),
      axis.text.x = element_text(size = 22, face = "bold.italic", 
                                 color = rep(my_colors[1:C], each = 3) ),
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 11) 
    ) +
    guides(color = "none")
  
  return(p)
}


lambda_plot <- function(data, n_sim, lwd, C, v, p,  
                        my_colors =  c("#6666FF", "#CC0066", "#009966", "#FF3300", "#660099"),
                        lty = "solid", v_just_ann = 0, h_just_ann = -1.0){
  labels <- get_lambda_labels(C = C, v = v, p = p)
  x_labels  <- labels$x_labels       # Returns bold(lambda[cvp])
  x_labb_hat<- labels$x_labb_hat     # Returns hat(lambda)[bold(cvp)]
  x_labb_real <- labels$x_labb_real    # Returns bold(lambda[cvp]^real)
  
  melted_data <- melt(data, id.vars = NULL)
  cl_ass <- rep(1:C, each = n_sim * p)
  data_for_plot <- cbind(melted_data, cl_ass)
  x <- seq(0.6, by = 1, length.out = C*p)
  xend <- seq(1.4, by = 1, length.out = C*p)
  x_ann <- seq(1, by = 1, length.out = C*p)
  y <- yend <- do.call(c, lapply(1:C, function(i) get(paste0("cat_probs", i))))
  p1 <- ggplot(data_for_plot, aes(x = variable, y = value, fill = as.factor(cl_ass), color = as.factor(cl_ass) ) ) +
    geom_boxplot(width = 0.5) 
  
  breaks_col <- NULL
  labels_col <- NULL
  values_col <- NULL
  for(i in 1:C){
    breaks_col <- c(breaks_col, as.character(i) )
    labels_col <- c(labels_col, as.character(i) )
    values_col <- c(values_col, setNames(my_colors[i], as.character(i) ))
  }
  ass_ <- rep(c(1:C), each = p)
  p1 <- p1 + purrr::map(1:(C*p), function(i) {
    geom_segment(aes(x = x[i], y = y[i], xend = xend[i], yend = yend[i], 
                     color = as.character(ass_[i])), linetype = lty, linewidth = lwd)
  } )
  
  p1 <- p1 + scale_color_manual(name = 'Cluster', breaks = breaks_col, 
                              labels = labels_col, values = alpha(values_col, 0.7))
  p1 <- p1 + scale_fill_manual(name = 'Cluster', breaks = breaks_col, 
                             labels = labels_col, values = alpha(values_col, 0.3))
  p1 <- p1 + labs(subtitle = "Lambda parameters",
                  title = "Categorical covariates",
         x = "Parameters",
         y = "Estimated Values", fill = "Cluster", color= "Cluster") +
    scale_x_discrete(labels = x_labb_hat) +
    annotate("text", x = x_ann, y = y, 
             label = x_labb_real, color = c(rep(my_colors[1:C], each = p)), 
             vjust = v_just_ann, hjust = h_just_ann, size = 5, parse = TRUE) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      plot.subtitle = element_text(size = 16, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold.italic"),
      legend.position = "right",
      legend.text = element_text(size = 18),
      legend.key.width = unit(3, "lines"),
      legend.title = element_text(face = "bold.italic", size = 18),
      axis.text.x = element_text(size = 22, face = "bold.italic", 
                                 color = rep(my_colors, each = 3) ),
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 11) 
    ) +
    guides(color = "none")
  
  return(p1)
}

inter_plot <- function(data, n_sim, lwd, C, vars,
                       my_colors = c("#6666FF", "#CC0066", "#009966", "#FF3300", "#660099"),
                       lty = "solid", v_just_ann = 0, h_just_ann = -1.0){
  u <- length(vars)
  labels <- get_inter_labels(C = C, dims = u)
  
  x_labels <- labels$x_labels     # Returns bold(gamma[cij])
  x_labb_hat <- labels$x_labb_hat    # Returns hat(gamma)[bold(cij)]
  x_labb_real <- labels$x_labb_real   # Returns bold(gamma[cij]^real)
  
  x <- seq(0.6, by = 1, length.out = C*(u*(u-1)/2) )
  xend <- seq(1.4, by = 1, length.out = C*(u*(u-1)/2) )
  y <- yend <- do.call(c, lapply(1:C_real, function(i) get(paste0("int", i))[upper.tri(get(paste0("int", i)))] ))
  
  melted_data <- melt(data, id.vars = NULL)
  cl_ass <- rep(1:C, each = n_sim * u*(u-1)/2)
  data_for_plot <- cbind(melted_data, cl_ass)
  
  p <- ggplot(data_for_plot, aes(x = variable, y = value, fill = as.factor(cl_ass), color = as.factor(cl_ass) ) ) +
    geom_boxplot(width = 0.5) 
  
  breaks_col <- NULL
  labels_col <- NULL
  values_col <- NULL
  for(i in 1:C){
    breaks_col <- c(breaks_col, as.character(i) )
    labels_col <- c(labels_col, as.character(i) )
    values_col <- c(values_col, setNames(my_colors[i], as.character(i) ))
  }
  ass_ <- rep(c(1:C), each = u*(u-1)/2)
  p <- p + purrr::map(1:(C*u*(u-1)/2), function(i) {
    geom_segment(aes(x = x[i], y = y[i], xend = xend[i], yend = yend[i], 
                     color = as.character(ass_[i])), linetype = lty, linewidth = lwd)
  } )
  
  p <- p + scale_color_manual(name = 'Cluster', breaks = breaks_col, 
                              labels = labels_col, values = alpha(values_col, 0.7))
  p <- p + scale_fill_manual(name = 'Cluster', breaks = breaks_col, 
                             labels = labels_col, values = alpha(values_col, 0.3))
  
  p <- p + labs(subtitle = "Interaction parameters",
                title = "Dichotomous dependent covariates",
                x = "Parameters",
                y = "Estimated Values", fill = "Cluster", color= "Cluster") +
    scale_x_discrete(labels = x_labb_hat) +
    annotate("text", x = seq(1, by = 1, length.out = C*(u*(u-1)/2) ), y = y, 
             label = x_labb_real, color = c(rep(my_colors[1:C], each = u*(u-1)/2)), 
             vjust = v_just_ann, hjust = h_just_ann, size = 5, parse = TRUE) +
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      plot.subtitle = element_text(size = 16, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold.italic"),
      legend.position = "right",
      legend.text = element_text(size = 18),
      legend.key.width = unit(3, "lines"),
      legend.title = element_text(face = "bold.italic", size = 18),
      axis.text.x = element_text(size = 22, face = "bold.italic", 
                                 color = rep(my_colors[1:C], each = 3) ),
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 11)
    ) +
    guides(color = "none")
  
  return(p)
}


thres_plot <- function(data, n_sim, lwd, C, vars,
                    my_colors = c("#6666FF", "#CC0066", "#009966", "#FF3300", "#660099"),
                    lty = "solid", v_just_ann = -0.2, h_just_ann = -1.8){
  labels <- get_thres_labels(C = C, variables = vars)
  u <- length(vars)
  x_labels <- labels$x_labels      # Returns bold(mu[c*x[var]])
  x_labb_hat <- labels$x_labb_hat    # Returns hat(mu)[bold(c*x[var])]
  x_labb_real <- labels$x_labb_real   # Returns bold(mu[c*x[var]]^real)
  
  x <- seq(0.6, by = 1, length.out = C*u)
  xend <- seq(1.4, by = 1, length.out = C*u)
  y <- yend <- do.call(c, lapply(1:C, function(i) get(paste0("thres", i))))
  
  melted_data <- melt(data, id.vars = NULL)
  cl_ass <- rep(1:C, each = n_sim * u)
  data_for_plot <- cbind(melted_data, cl_ass)
  
  p <- ggplot(data_for_plot, aes(x = variable, y = value, fill = as.factor(cl_ass), color = as.factor(cl_ass) ) ) +
    geom_boxplot(width = 0.5) 
  
  breaks_col <- NULL
  labels_col <- NULL
  values_col <- NULL
  for(i in 1:C){
    breaks_col <- c(breaks_col, as.character(i) )
    labels_col <- c(labels_col, as.character(i) )
    values_col <- c(values_col, setNames(my_colors[i], as.character(i) ))
  }
  ass_ <- rep(c(1:C), each = u)
  p <- p + purrr::map(1:(C*u), function(i) {
    geom_segment(aes(x = x[i], y = y[i], xend = xend[i], yend = yend[i], 
                     color = as.character(ass_[i])), linetype = lty, linewidth = lwd)
  } )
  
  p <- p + scale_color_manual(name = 'Cluster', breaks = breaks_col, 
                              labels = labels_col, values = alpha(values_col, 0.7))
  p <- p + scale_fill_manual(name = 'Cluster', breaks = breaks_col, 
                             labels = labels_col, values = alpha(values_col, 0.3))
  p <- p +
    labs(subtitle = "Threshold parameters",
         title = "Dichotomous dependent covariates",
         x = "Parameters",
         y = "Estimated Values", fill = "Cluster", color= "Cluster") +
    scale_x_discrete(labels = x_labb_hat) +
    annotate("text", x = seq(1, by = 1, length.out = C*u), y = y, 
             label = x_labb_real, color = c(rep(my_colors[1:C], each = u)), 
             vjust = v_just_ann, hjust = h_just_ann, size = 5, parse = TRUE) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      plot.subtitle = element_text(size = 16, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold.italic"),
      legend.position = "right",
      legend.text = element_text(size = 18),
      legend.key.width = unit(3, "lines"),
      legend.title = element_text(face = "bold.italic", size = 18),
      axis.text.x = element_text(size = 22, face = "bold.italic", 
                                 color = rep(my_colors[1:C], each = u) ),
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 11) 
    ) +
    guides(color = "none")
  
  return(p)
}

#------------------------------Support functions--------------------------------
get_beta_labels <- function(variables) {
  n_vars <- length(variables)
  x_labels <- vector("list", n_vars)
  x_labb_hat <- vector("list", n_vars)
  x_labb_real <- vector("list", n_vars)
  
  for (var_index in 1:n_vars) {
      cluster_char <- as.character("c")
      idx <- var_index
      var_name <- variables[var_index]
      x_labels[[idx]] <- bquote(bold(beta[.(cluster_char) * .(var_name)]))
      x_labb_hat[[idx]] <- bquote(hat(beta)[bold(.(cluster_char) * .(var_name))])
      x_labb_real[[idx]] <- bquote(bold(beta[.(cluster_char) * .(var_name)]^real))
  }
  
  return(list(x_labels = x_labels, x_labb_hat = x_labb_hat, x_labb_real = x_labb_real))
}


get_lambda_labels <- function(C, v, p) {
  x_labels <- vector("list", C * p)
  x_labb_hat <- vector("list", C * p)
  x_labb_real <- vector("list", C * p)
  
  for (cluster in 1:C) {
    for (param_index in 1:p) {
      idx <- (cluster - 1) * p + param_index
      x_labels[[idx]] <- bquote(bold(lambda[.(cluster * 100 + v * 10 + param_index)]))
      x_labb_hat[[idx]] <- bquote(hat(lambda)[bold(.(cluster * 100 + v * 10 + param_index))])
      x_labb_real[[idx]] <- bquote(bold(lambda[.(cluster * 100 + v * 10 + param_index)]^real))
    }
  }
  return(list(x_labels = x_labels, x_labb_hat = x_labb_hat, x_labb_real = x_labb_real))
}

get_mu_labels <- function(C, variables) {
  n_vars <- length(variables)
  x_labels <- vector("list", C * n_vars)
  x_labb_hat <- vector("list", C * n_vars)
  x_labb_real <- vector("list", C * n_vars)
  
  for (cluster in 1:C) {
    for (var_index in 1:n_vars) {
      cluster_char <- as.character(cluster)
      idx <- (cluster - 1) * n_vars + var_index
      var_name <- variables[var_index]
      x_labels[[idx]] <- bquote(bold(mu[.(cluster_char) * x[.(var_name)]]))
      x_labb_hat[[idx]] <- bquote(hat(mu)[bold(.(cluster_char) * x[.(var_name)])])
      x_labb_real[[idx]] <- bquote(bold(mu[.(cluster_char) * x[.(var_name)]]^real))
    }
  }
  return(list(x_labels = x_labels, x_labb_hat = x_labb_hat, x_labb_real = x_labb_real))
}

get_sigma_labels <- function(C, dims) {
  x_labels <- vector("list", C * (dims * (dims + 1) / 2))
  x_labb_hat <- vector("list", C * (dims * (dims + 1) / 2))
  x_labb_real <- vector("list", C * (dims * (dims + 1) / 2))
  
  idx <- 1
  for (cluster in 1:C) {
    for (i in 1:dims) {
      for (j in i:dims) {
        i_char <- as.character(i)
        j_char <- as.character(j)
        cluster_char <- as.character(cluster)

        x_labels[[idx]] <- bquote(bold(Sigma[.(cluster_char) * .(i_char) * .(j_char)]))
        x_labb_hat[[idx]] <- bquote(hat(Sigma)[bold(.(cluster_char) * .(i_char) * .(j_char))])
        x_labb_real[[idx]] <- bquote(bold(Sigma[.(cluster_char) * .(i_char) * .(j_char)]^real))
        
        # Increment index
        idx <- idx + 1
      }
    }
  }
  
  return(list(x_labels = x_labels, x_labb_hat = x_labb_hat, x_labb_real = x_labb_real))
}

get_inter_labels <- function(C, dims) {
  x_labels <- vector("list", C * ((dims * (dims - 1)) / 2))
  x_labb_hat <- vector("list", C * ((dims * (dims - 1)) / 2))
  x_labb_real <- vector("list", C * ((dims * (dims - 1)) / 2))
  
  idx <- 1
  idx <- 1
  for (cluster in 1:C) {
    for (i in 1:(dims - 1)) {
      for (j in (i + 1):dims) {
        i_char <- as.character(i)
        j_char <- as.character(j)
        cluster_char <- as.character(cluster)
        
        x_labels[[idx]] <- bquote(bold(gamma[.(cluster_char) * .(i_char) * .(j_char)]))
        x_labb_hat[[idx]] <- bquote(hat(gamma)[bold(.(cluster_char) * .(i_char) * .(j_char))])
        x_labb_real[[idx]] <- bquote(bold(gamma[.(cluster_char) * .(i_char) * .(j_char)]^real))
        
        idx <- idx + 1
      }
    }
  }
  return(list(x_labels = x_labels, x_labb_hat = x_labb_hat, x_labb_real = x_labb_real))
}

get_thres_labels <- function(C, variables) {
  n_vars <- length(variables)
  x_labels <- vector("list", C * n_vars)
  x_labb_hat <- vector("list", C * n_vars)
  x_labb_real <- vector("list", C * n_vars)
  
  for (cluster in 1:C) {
    for (var_index in 1:n_vars) {
      idx <- (cluster - 1) * n_vars + var_index
      var_name <- variables[var_index]
      clust <- as.character(cluster)
      x_labels[[idx]] <- bquote(bold(nu[.(clust) * .(var_name)]))
      x_labb_hat[[idx]] <- bquote(hat(nu)[bold(.(clust) * .(var_name))])
      x_labb_real[[idx]] <- bquote(bold(nu[.(clust) * .(var_name)]^real))
    }
  }
  return(list(x_labels = x_labels, x_labb_hat = x_labb_hat, x_labb_real = x_labb_real))
}
  
upper_diag_strings <- function(N) {
  
  mat <- matrix(1:(N*N), nrow=N, ncol=N)
  indices <- which(row(mat) <= col(mat), arr.ind = TRUE)
  strings <- paste(indices[, 1], indices[, 2], sep = "")
  
  return(strings)
}

get_comp_models_beta <- function(fit_glm, fit_glmer){
  beta_glm <- data.frame(matrix(0,nrow = K, ncol = length(coef(fit_glm[[1]]))-1 ))
  beta_glmer <- data.frame(matrix(0,nrow = K, ncol = length(fixef(fit_glmer[[1]]))-1 ))
  for(k in 1:K){
    beta_glm[k,] <- coef(fit_glm[[k]])[2:(length(coef(fit_glm[[1]])))]
    beta_glmer[k,] <- fixef(fit_glmer[[k]])[2:(length(fixef(fit_glmer[[1]])))]
  }
  comp_data <- rbind(beta_glm, beta_glmer)
  return(comp_data)
}


split_dataset_by_columns <- function(data, N, column_names) {
  total_columns <- ncol(data)
  
  if (length(column_names) != N) {
    stop("The length of column_names must be equal to N.")
  }
  num_datasets <- ceiling(total_columns / N)
  subdatasets <- list()
  for (i in 1:num_datasets) {
    start_col <- ((i - 1) * N) + 1
    end_col <- min(i * N, total_columns)
    subdata <- data[, start_col:end_col]
    colnames(subdata) <- column_names[1:ncol(subdata)]
    subdatasets[[i]] <- subdata
  }
  return(subdatasets)
}



