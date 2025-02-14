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

K <- 100
sim_data <- list()
for(i in 1:K){
  sim_data[[i]] <- readRDS( paste0("sim_data/mydata", i, ".rds") )
}
test_data <- readRDS( "sim_data/testdata.rds" )

# Set real number of clusters
C_real <- 2

#----------------------------------FLEXMIX--------------------------------------
library(flexmix)

yy <- cbind(success = sim_data[[1]]$y, failure = 1 - sim_data[[1]]$y)
fit_fm <- flexmix(yy ~ x1 + x2 + factor(V1) + factor(D1) + factor(D2) + factor(D3), 
                       data = sim_data[[1]], k = 2, model = FLXMRglm(family = "binomial"))


fit_fm <- list()
set.seed(77)
for(i in 1:K){
  yy <- cbind(success = sim_data[[i]]$y, failure = 1 - sim_data[[i]]$y)
  init_clusters <- kmeans(sim_data[[i]][, c("x1", "x2")], centers = 2)$cluster
  # MULTIPLE RANDOM INIT
  # fit_fm <- flexmix(yy ~ x1 + x2 + factor(V1) + factor(D1) + factor(D2) + factor(D3), 
  #                    data = sim_data[[i]], k = 2, model = FLXMRglm(family = "binomial"),
  #                   control = list(nrep = 5))
  # K-MEANS INIT
  fit_fm[[i]] <- flexmix(yy ~ x1 + x2 + factor(V1) + factor(D1) + factor(D2) + factor(D3), 
                    data = sim_data[[i]], k = 2, model = FLXMRglm(family = "binomial"),
                    control = list(nrep = 10), cluster = init_clusters)
  print(i)
}

ari_flexmix <- vector(mode = "numeric", length = K)
miss_flexmix <- vector(mode = "numeric", length = K)
acc_train <- vector(mode = "numeric", length = K)
acc_test <- vector(mode = "numeric", length = K)
ece <- vector(mode = "numeric", length = K)
brier <- vector(mode = "numeric", length = K)

for(i in 1:K){
  fitted_clusters <- fit_fm[[i]]@cluster
  true_clusters <- sim_data[[i]]$latent
  tab_obs <- table(fitted_clusters, true_clusters)
  perm <- solve_LSAP(tab_obs, maximum = TRUE) 
  
  matched_confusion <- tab_obs[cbind(1:nrow(tab_obs), perm)]
  misclassified <- sum(tab_obs) - sum(matched_confusion) 
  miss_flexmix[i] <- misclassified / sum(tab_obs)
  ari_flexmix[i] <- adjustedRandIndex(true_clusters, fitted_clusters)
  print(i)
  
  predicted_probs_train <- predict(fit_fm[[i]], aggregate = T)[[1]][,1]
  roc_curve <- roc(sim_data[[i]][,"y"], predicted_probs_train, levels = c(0, 1),
                   direction = "<")
  optimal.threshold <- coords(roc_curve, "best", best.method = "closest.topleft")
  predicted_classes_train <- ifelse(predicted_probs_train > as.numeric(optimal.threshold[1]), 1, 0)
  acc_train[i] <- mean(predicted_classes_train == sim_data[[i]][,"y"])
  
  predicted_probs_test <- predict(fit_fm[[i]], newdata = test_data, aggregate = T)[[1]][,1]
  predicted_classes_test <- ifelse(predicted_probs_test > as.numeric(optimal.threshold[1]), 1, 0)
  acc_test[i] <- mean(predicted_classes_test == test_data[,"y"])
  
  ece[i] <- getECE(test_data$y, predicted_probs_test)
  brier[i] <- get_Brier_score(test_data$y, predicted_probs_test)$brier
}

mean(ari_flexmix)
sd(ari_flexmix)

mean(miss_flexmix)
sd(miss_flexmix)

mean(acc_train)
sd(acc_train)

mean(acc_test)
sd(acc_test)

mean(ece)
sd(ece)

mean(brier)
sd(brier)

#----------------------------------MIX_TOOLS------------------------------------
library(mixtools)

design_matrix <- model.matrix(~ x1 + x2 + factor(V1) + factor(D1) + factor(D2) 
                              + factor(D3) -1 , data = sim_data[[1]])

fit <- logisregmixEM(y = sim_data[[1]]$y, x = design_matrix, k = 2,
                     addintercept = TRUE)

# WORKS ONLY FOR a SUBSET
design_matrix <- model.matrix(~ x1 + x2 
                              -1 , data = sim_data[[1]])
fit <- logisregmixEM(y = sim_data[[1]]$y, x = design_matrix, k = 2,
                     addintercept = TRUE)

summary(fit)
print(fit$lambda)

print(fit$beta)
head(fit$posterior)
clusters <- apply(fit$posterior, 1, which.max)
