### Import and set directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  path <- getwd()
  print(getwd())
} else {
  message("The rstudioapi package is not available.")
}
source(paste0( path,"/visual_functions.R") )
source( paste0(dirname(dirname(path)), "/required_packages.R") )
source( paste0(dirname(dirname(path)), "/Functions/EM-algorithm.R") )
source( paste0(dirname(dirname(path)), "/Functions/support_functions.R") )
source( paste0(dirname(dirname(path)), "/Functions/initialization_functions.R") )
source( paste0(dirname(dirname(path)), "/Functions/E-step.R") )
source( paste0(dirname(dirname(path)), "/Functions/M-step.R") )
source( paste0(dirname(dirname(path)), "/Functions/log_likelihood.R") )
source( paste0(dirname(dirname(path)), "/Functions/BIC_calculation.R") )
source( paste0(dirname(dirname(path)), "/Functions/Predictions.R") )

# Load simulated datasets and model fitting
fit_list_full = readRDS("fitted_models/fit_list_full")
fit_list_partial = readRDS("fitted_models/fit_list_partial")

fit_glm = readRDS("fitted_models/fit_list_glm")
fit_glmer = readRDS("fitted_models/fit_list_glmer")

K <- 100
sim_data <- list()
for(i in 1:K){
  sim_data[[i]] <- readRDS( paste0("sim_data/mydata", i, ".rds") )
}
test_data <- readRDS( "sim_data/testdata.rds" )

# Set real number of clusters
C_real <- 2

# Visualization of clusters on continuous variables space (e.g., simulated data 2)
p_real_clust <- cluster_plot_byclust(sim_data[[10]], sim_data[[10]]$latent,
                                     add_to_title = " (Real clusters)")
p_fitted_clust <- cluster_plot(sim_data[[10]], fit_list_full[[C_real]][[10]],
                               add_to_title = " (Fitted clusters)")
grid.arrange(p_real_clust, p_fitted_clust, nrow = 1, ncol = 2)

#-------------------------------BIC comparison----------------------------------
bic_full <- list()
full_names <- paste0("ML-CWMd \n(C = ", 1:3, ")")
for(k in 1:3){
  bic_full[[k]] <- unlist(lapply(fit_list_full[[k]], function(x) x$bic))
  
} 
names(bic_full) <- full_names

bic_partial <- list()
partial_names <- paste0("ML-CWMd \nwith no D \n(C = ", 1:3, ")")
for(k in 1:3){
  bic_partial[[k]] <- unlist(lapply(fit_list_partial[[k]], function(x) x$bic))
} 
names(bic_partial) <- partial_names

bic_list <- c(bic_full, bic_partial)
p_bic <- BIC_comparison_plot(bic_list = bic_list, n_sim = K, 
                lwd = 1, lty = "dashed", 
                v_just_ann = -0.2, h_just_ann = -0.7,
                my_colors = c("#6666FF", "#CC0066", "#009966") )

min_indices <- sapply(1:length(bic_full[[1]]), function(i) {
  # Get the vector of the current index across all vectors
  current_values <- sapply(bic_full, function(vec) vec[i])
  # Find the index of the minimum value
  which.min(current_values)  # Returns the index of the minimum value
})
#------------------------------Fitted vs Real clusters--------------------------

# Get percentage of misclassified observations
miss_full <- calculate_misclassification(fit_list = fit_list_full, 
             sim_data_list = sim_data, C_real = C_real, n_sim = K)
miss_partial <- calculate_misclassification(fit_list = fit_list_partial, 
                sim_data_list = sim_data, C_real = C_real, n_sim = K) 

# Get ARI
ari_full <- calculate_ARI(fit_list = fit_list_full, 
                          sim_data_list = sim_data, C_real = C_real, n_sim = K)
ari_partial <- calculate_ARI(fit_list = fit_list_partial, 
                             sim_data_list = sim_data, C_real = C_real, n_sim = K)

p1 <- MLCWMd_comparison(full = miss_full, partial = miss_partial, 
                        x_breaks = seq(0.0, 1, 0.03), title = "Percentage of Misclassified Observations", 
                        y_title = "Percentage (%)", col_x_axis = F, show.outliers = F)

p2 <- MLCWMd_comparison(full = ari_full, partial = ari_partial, 
                        x_breaks = seq(0.0, 1, 0.05), title = "Adjusted Rand Index Comparison", 
                        y_title = "ARI", col_x_axis = F, show.outliers = F)
grid.arrange(p1, p2, ncol = 2, nrow = 1)

# Get Calibration
ece_full <- get_MLCWMd_calib(fit_list = fit_list_full, 
                             test_data = test_data, C_real = C_real, n_sim = K, flag = F)
ece_partial <- get_MLCWMd_calib(fit_list = fit_list_partial, 
                                test_data = test_data, C_real = C_real, n_sim = K, flag = F)

ece_glm <- get_glm_calib(fit_list = fit_glm, 
                         test_data = test_data, n_sim = K, flag = F)
ece_glmer <- get_glmer_calib(fit_list = fit_glmer, 
                             test_data = test_data, n_sim = K, flag = F)

p3 <- MLCWMd_comparison(full = ece_full["ECE"]$ECE, partial = ece_partial["ECE"]$ECE,
                        glm_res = ece_glm["ECE"]$ECE, glmer_res = ece_glmer["ECE"]$ECE,
                        x_breaks = seq(0.0, 1, 0.03), title = "ECE on Test Data", 
                        y_title = "ECE", col_x_axis = F, show.outliers = F)
p4 <- MLCWMd_comparison(full = ece_full["MCE"]$MCE, partial = ece_partial["MCE"]$MCE,
                        glm_res = ece_glm["MCE"]$MCE, glmer_res = ece_glmer["MCE"]$MCE,
                        x_breaks = seq(0.0, 1, 0.01), title = "MCE Index on Test Data", 
                        y_title = "MCE", col_x_axis = T, show.outliers = F)
p5 <- MLCWMd_comparison(full = ece_full["Brier"]$Brier, partial = ece_partial["Brier"]$Brier,
                        glm_res = ece_glm["Brier"]$Brier, glmer_res = ece_glmer["Brier"]$Brier,
                        x_breaks = seq(0.0, 1, 0.03), title = "Brier Score on Test Data", 
                        y_title = "Brier Score", col_x_axis = F, show.outliers = F)
grid.arrange(p3, p4, p5, ncol = 2, nrow = 2)

#-------------------------Train&Test Prediction Accuracy------------------------
pred_acc_full <- get_MLCWMd_pred_accuracy(fit_list = fit_list_full, 
          sim_data_list = sim_data, test_data = test_data, C_real = C_real, n_sim = K)

pred_acc_partial <- get_MLCWMd_pred_accuracy(fit_list = fit_list_partial, 
                          sim_data_list = sim_data, test_data = test_data, 
                                      C_real = C_real, n_sim = K)

pred_acc_glm <- get_glm_pred_acc(fit_list = fit_glm, 
                        sim_data_list = sim_data, test_data = test_data, n_sim = K)

pred_acc_glmer <- get_glm_pred_acc(fit_list = fit_glmer, 
                      sim_data_list = sim_data, test_data = test_data, n_sim = K)

p_train <- MLCWMd_comparison(full = pred_acc_full[["train_acc"]], 
                             partial = pred_acc_partial[["train_acc"]],
                             glm_res = pred_acc_glm[["train_acc"]], 
                             glmer_res = pred_acc_glmer[["train_acc"]],
                             x_breaks = seq(0.60, 0.85, 0.03), 
                             title = "Predicted Responses Accuracy on Training Data",
                             y_title = "Accuracy", col_x_axis = F)

p_test <- MLCWMd_comparison(full = pred_acc_full[["test_acc"]], 
                            partial = pred_acc_partial[["test_acc"]],
                            glm_res = pred_acc_glm[["test_acc"]], 
                            glmer_res = pred_acc_glmer[["test_acc"]],
                            x_breaks = seq(0.55, 0.85, 0.03), 
                            title = "Predicted Responses Accuracy on Test Data",
                            y_title = "Accuracy", col_x_axis = F)

p_train <- MLCWMd_comparison(full = pred_acc_full[["train_acc"]], 
                             partial = pred_acc_partial[["train_acc"]],
                             x_breaks = seq(0.60, 0.85, 0.01), 
                             title = "Predicted Responses Accuracy on Training Data",
                             y_title = "Accuracy", col_x_axis = T)

p_test <- MLCWMd_comparison(full = pred_acc_full[["test_acc"]], 
                            partial = pred_acc_partial[["test_acc"]],
                            x_breaks = seq(0.60, 0.85, 0.02), 
                            title = "Predicted Responses Accuracy on Test Data",
                            y_title = "Accuracy", col_x_axis = T)
grid.arrange(p_train, p_test, nrow = 1, ncol = 2)
#grid.arrange(p_train, p_test, p3, p6, ncol = 2, nrow = 2)


#--------------------------Visualize parameters recovery------------------------
load_real_params()

#--------------------Multilevel Logistic Regression Part------------------------
### BETA parameters
beta_data <- extract_beta_from_sim_data(fit_list = fit_list_full, 
             sim_data_list = sim_data, C = 2)
comp_data <- get_comp_models_beta(fit_glm, fit_glmer)

p_beta <- beta_plot(data = beta_data, comp_data = comp_data, n_sim = dim(beta_data)[1], 
                    lwd = 0.6, C = 2, vars = c("1", "2", "3", "4", "5",
                                               "6", "7"), lty = "dashed", 
                    v_just_ann = -0.2, h_just_ann = -0.7, rm.outliers = FALSE)

beta_data <- extract_beta_from_sim_data(fit_list = fit_list_partial, 
                                        sim_data_list = sim_data, C = 2)
p_beta_partial <- beta_plot(data = beta_data, comp_data = comp_data, n_sim = dim(beta_data)[1], 
                            lwd = 0.6, C = 2, vars = c("1", "2", "3", "4", "5",
                                                       "6", "7"), lty = "dashed", 
                            v_just_ann = -0.2, h_just_ann = -0.7, rm.outliers = FALSE)
p_beta

#------------Continuous covariates parameters (mu and sigma)--------------------
### MU parameters
mu_data <- extract_mu_from_sim_data(fit_list = fit_list_full, 
                sim_data_list = sim_data, C = 2)
p_mu <- mu_plot(data = mu_data, n_sim = dim(mu_data)[1], 
      lwd = 1, C = 2, vars = c("1", "2"), lty = "dashed", 
      v_just_ann = -0.05, h_just_ann = -1.7)

mu_data <- extract_mu_from_sim_data(fit_list = fit_list_partial, 
           sim_data_list = sim_data, C = 2)
p_mu_partial <- mu_plot(data = mu_data, n_sim = dim(mu_data)[1],
      lwd = 1, C = 2, vars = c("1", "2"), lty = "dashed", 
      v_just_ann = -0.05, h_just_ann = -1.7)

p_mu 
p_mu_partial

### SIGMA parameters
# Full
sigma_data <- extract_sigma_from_sim_data(fit_list = fit_list_full, 
              sim_data_list = sim_data, C = 2)

p_sigma <- sigma_plot(data = sigma_data, n_sim = dim(sigma_data)[1], lwd = 1, 
      C = 2, vars = c("1", "2"), lty = "dashed", v_just_ann = -0.2, h_just_ann = -1.15)

# Partial
sigma_data <- extract_sigma_from_sim_data(fit_list = fit_list_partial, 
              sim_data_list = sim_data, C = 2)

p_sigma_partial <- sigma_plot(data = sigma_data, n_sim = dim(sigma_data)[1], lwd = 1, 
      C = 2, vars = c("1", "2"), lty = "dashed", v_just_ann = -0.2, h_just_ann = -1.15)
p_sigma
p_sigma_partial

#-----------------Categorical covariates parameters (lambda)--------------------
### LAMBDA parameters
# Full
lam_data <- extract_lambda_from_sim_data(fit_list = fit_list_full, 
            sim_data_list = sim_data, C = 2, v = 1)
p_lam <- lambda_plot(data = lam_data, n_sim = dim(lam_data)[1], 
      lwd = 1, C = 2, v = 1, p = 3, lty = "dashed", v_just_ann = -0.2, h_just_ann = -1.15)

# Partial
lam_data <- extract_lambda_from_sim_data(fit_list = fit_list_partial, 
                                         sim_data_list = sim_data, C = 2, v = 1)
p_lam_partial <- lambda_plot(data = lam_data, n_sim = dim(lam_data)[1], 
      lwd = 1, C = 2, v = 1, p = 3, lty = "dashed", v_just_ann = -0.2, h_just_ann = -1.15)

p_lam
p_lam_partial

#------------------Dichotomous dependent covariates-----------------------------
### INTERACTIONS parameters
int_data <- extract_inter_from_sim_data(fit_list = fit_list_full, 
                                        sim_data_list = sim_data, C = 2)

p_inter <- inter_plot(data = int_data, n_sim = dim(int_data)[1], lwd = 1, 
                 C = 2, vars = c("1", "2", "3"),
                 lty = "dashed", v_just_ann = -0.2, h_just_ann = -1.2)
p_inter

### THRESHOLD parameters
thres_data <- extract_thres_from_sim_data(fit_list = fit_list_full, 
              sim_data_list = sim_data, C = 2)
p_thres <- thres_plot(data = thres_data, n_sim = dim(thres_data)[1], 
              lwd = 1, C = 2, vars = c("1", "2", "3"), lty = "dashed",
               v_just_ann = -0.2, h_just_ann = -1.2)
p_thres
