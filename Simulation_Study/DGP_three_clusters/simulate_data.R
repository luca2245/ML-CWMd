if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  path <- getwd()
  print(getwd())
} else {
  message("The rstudioapi package is not available.")
}

# Load necessary packages
source( paste0(dirname(dirname(path)), "/required_packages.R") )
source( paste0(dirname(dirname(path)), "/Functions/support_functions.R") )


R=100
dir.create("sim_data")

for (i in 1:R){
  
  set.seed(777*i)
  
  filename_data <- paste("sim_data/mydata", i, ".rds", sep = "")
  
  # Set DPG parameters
  n_obs <- 1500
  n_groups <- 10 # number of groups
  n_per_group <- 150 # number of observations per group
  group_ids <- rep(1:n_groups, each = n_per_group)
  w <- c(0.20,0.30,0.50) # proportion of observations per cluster
  
  # Random Intercept
  group_intercepts_1 <- c(-0.11, -0.38, -0.26, -2.83, -0.79, -0.76, 0.71, 0.93, -2.15, -1.46)
  group_intercepts_2 <- c(1.16, 0.52, -1.05, -0.10, 0.61, 1.62, -0.42, 0.41, 0.05, -1.45)
  group_intercepts_3 <- c(-1.18, -1.05, -0.20, -0.95, -0.84, 1.85, -0.55, 0.21, -1.83, -1.49)
  
  # MU and SIGMA Parameters
  cov.1 <- matrix(c(1,0.5,0.5,1),nrow = 2)
  mean.1 <- c(2, 0)
  
  cov.2 <- matrix(c(1,0.5,0.5,1),nrow = 2)
  mean.2 <- c(4.5, 3.5)
  
  cov.3 <- matrix(c(1,0.5,0.5,1),nrow = 2)
  mean.3 <- c(4, 0)
  
  # LAMBDA Parameters
  cat_probs1 <- c(0.7, 0.2, 0.1)
  cat_probs2 <- c(0.15, 0.65, 0.2)
  cat_probs3 <- c(0.3, 0.5, 0.2)
  
  #INTERACTION AND THRESHOLD Parameters
  N <- 3
  
  int1 <- matrix(c(0, 1.8, -2.1, 1.8, 0, -1.4, -2.1, -1.4, 0), nrow = 3, ncol = 3)
  thres1 <- c(0, 0, 0)
  
  int2 <- matrix(c(0, -1.5, 2.8, -1.5,  0, 0, 2.8,  0, 0), nrow = 3, ncol = 3)
  thres2 <- c(0, 0, 0)
  
  int3 <- matrix(c(0, -4.2, 1.1, -4.2, 0,  1.3, 1.1,  1.3,  0), nrow = 3, ncol = 3)
  thres3 <- c(0, 0, 0)
  
  # Fixed effects
  beta1 <- c(-1,  0.5, -1.3,  0.2,  1.9,  2.7, -2.3) 
  beta2 <- c(-0.2,  0.7, -0.4,  0.4, -2.5, -0.6,  0.6)
  beta3 <- c(-0.5, -0.5, 1.3, -0.6, -2.0,  1.8,  3.2)

  # Simulate data
  data <- as.data.frame(rbind(mvrnorm(n = w[1]*n_obs, mu = mean.1, Sigma = cov.1), mvrnorm(n = w[2]*n_obs, mu = mean.2, Sigma = cov.2), mvrnorm(n = w[3]*n_obs, mu = mean.3, Sigma = cov.3)))
  data <- cbind(rep(1:n_groups, each = n_per_group)[sample(n_obs)],
                c(rep(1, length.out = w[1]*n_obs),rep(2, length.out = w[2]*n_obs),
                  rep(3, length.out = w[3]*n_obs)), data)
  colnames(data) <- c('level','latent','x1','x2')
  
  result1 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = cat_probs1),
                   rmultinom(w[2]*n_obs, size = 1, prob = cat_probs2),
                   rmultinom(w[3]*n_obs, size = 1, prob = cat_probs3))
  
  V <- data.frame( v1 = data.frame(t(result1)) )
  colnames(V) <- c("A","B","C")
  data <- cbind(data,V)
  
  dich.1 <- as.data.frame(IsingSampler(w[1]*n_obs,int1,thres1,beta = 1))
  dich.2 <- as.data.frame(IsingSampler(w[2]*n_obs,int2,thres2, beta = 1))
  dich.3 <- as.data.frame(IsingSampler(w[3]*n_obs,int3,thres3, beta = 1))
  dich <- rbind(dich.1,dich.2,dich.3)
  colnames(dich) <- c("D1","D2","D3")
  data <- cbind(data,dich)

  data <- data %>%
    dplyr::mutate(lin_pred = 0) %>%
    dplyr::mutate(lin_pred = if_else(latent == 1,
                                     group_intercepts_1[level] + beta1[1]*x1 + beta1[2]*x2 + beta1[3]*B + beta1[4]*C + 
                                       beta1[5]*D1 + beta1[6]*D2 + beta1[7]*D3,lin_pred)) %>%
    dplyr::mutate(lin_pred = if_else(latent == 2, 
                                     group_intercepts_2[level] + beta2[1]*x1 + beta2[2]*x2 + beta2[3]*B + beta2[4]*C + 
                                       beta2[5]*D1 + beta2[6]*D2 + beta2[7]*D3, lin_pred)) %>%
    dplyr::mutate(lin_pred = if_else(latent == 3,
                                     group_intercepts_3[level] +  beta3[1]*x1 + beta3[2]*x2 + beta3[3]*B + beta3[4]*C + 
                                     beta3[5]*D1 + beta3[6]*D2 + beta3[7]*D3, lin_pred ))


  y <- vector("numeric",dim(data)[1])
  for(i in 1:dim(data)[1]){
    y[i] <- rbinom(1, size = 1, prob = plogis(data$lin_pred[i]))
  }
  data <- cbind(y,data)
  data <- data[, !names(data) %in% 'lin_pred']
  
  data <- data[sample(nrow(data)),]
  data = onehot_to_cat(data, list(c("A", "B", "C")), "V1")
  
  check.1 <- sum(data$y[which(data$latent == 1)])/(w[1]*n_obs)
  check.2 <- sum(data$y[which(data$latent == 2)])/(w[2]*n_obs)
  check.3 <- sum(data$y[which(data$latent == 3)])/(w[3]*n_obs)
  
  saveRDS(data, file=filename_data)

}

#------------------------------Get Test Data------------------------------------
set.seed(22)
filename_data <- paste("sim_data/testdata", ".rds", sep = "")

# Set DPG parameters
n_obs <- 1500
n_groups <- 10 # number of groups
n_per_group <- 150 # number of observations per group
group_ids <- rep(1:n_groups, each = n_per_group)
w <- c(0.20,0.30,0.50) # proportion of observations per cluster

# Random Intercept
group_intercepts_1 <- c(-0.11, -0.38, -0.26, -2.83, -0.79, -0.76, 0.71, 0.93, -2.15, -1.46)
group_intercepts_2 <- c(1.16, 0.52, -1.05, -0.10, 0.61, 1.62, -0.42, 0.41, 0.05, -1.45)
group_intercepts_3 <- c(-1.18, -1.05, -0.20, -0.95, -0.84, 1.85, -0.55, 0.21, -1.83, -1.49)

# MU and SIGMA Parameters
cov.1 <- matrix(c(1,0.5,0.5,1),nrow = 2)
mean.1 <- c(2, 0)

cov.2 <- matrix(c(1,0.5,0.5,1),nrow = 2)
mean.2 <- c(4.5, 3.5)

cov.3 <- matrix(c(1,0.5,0.5,1),nrow = 2)
mean.3 <- c(4, 0)

# LAMBDA Parameters
cat_probs1 <- c(0.7, 0.2, 0.1)
cat_probs2 <- c(0.15, 0.65, 0.2)
cat_probs3 <- c(0.3, 0.5, 0.2)

#INTERACTION AND THRESHOLD Parameters
N <- 3

int1 <- matrix(c(0, 1.8, -2.1, 1.8, 0, -1.4, -2.1, -1.4, 0), nrow = 3, ncol = 3)
thres1 <- c(0, 0, 0)

int2 <- matrix(c(0, -1.5, 2.8, -1.5,  0, 0, 2.8,  0, 0), nrow = 3, ncol = 3)
thres2 <- c(0, 0, 0)

int3 <- matrix(c(0, -4.2, 1.1, -4.2, 0,  1.3, 1.1,  1.3,  0), nrow = 3, ncol = 3)
thres3 <- c(0, 0, 0)

# Fixed effects
beta1 <- c(-1,  0.5, -1.3,  0.2,  1.9,  2.7, -2.3) 
beta2 <- c(-0.2,  0.7, -0.4,  0.4, -2.5, -0.6,  0.6)
beta3 <- c(-0.5, -0.5, 1.3, -0.6, -2.0,  1.8,  3.2)

# Simulate data
data <- as.data.frame(rbind(mvrnorm(n = w[1]*n_obs, mu = mean.1, Sigma = cov.1), mvrnorm(n = w[2]*n_obs, mu = mean.2, Sigma = cov.2), mvrnorm(n = w[3]*n_obs, mu = mean.3, Sigma = cov.3)))
data <- cbind(rep(1:n_groups, each = n_per_group)[sample(n_obs)],
              c(rep(1, length.out = w[1]*n_obs),rep(2, length.out = w[2]*n_obs),
                rep(3, length.out = w[3]*n_obs)), data)
colnames(data) <- c('level','latent','x1','x2')

result1 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = cat_probs1),
                 rmultinom(w[2]*n_obs, size = 1, prob = cat_probs2),
                 rmultinom(w[3]*n_obs, size = 1, prob = cat_probs3))

V <- data.frame( v1 = data.frame(t(result1)) )
colnames(V) <- c("A","B","C")
data <- cbind(data,V)

dich.1 <- as.data.frame(IsingSampler(w[1]*n_obs,int1,thres1,beta = 1))
dich.2 <- as.data.frame(IsingSampler(w[2]*n_obs,int2,thres2, beta = 1))
dich.3 <- as.data.frame(IsingSampler(w[3]*n_obs,int3,thres3, beta = 1))
dich <- rbind(dich.1,dich.2,dich.3)
colnames(dich) <- c("D1","D2","D3")
data <- cbind(data,dich)

data <- data %>%
  dplyr::mutate(lin_pred = 0) %>%
  dplyr::mutate(lin_pred = if_else(latent == 1,
                                   group_intercepts_1[level] + beta1[1]*x1 + beta1[2]*x2 + beta1[3]*B + beta1[4]*C + 
                                     beta1[5]*D1 + beta1[6]*D2 + beta1[7]*D3,lin_pred)) %>%
  dplyr::mutate(lin_pred = if_else(latent == 2, 
                                   group_intercepts_2[level] + beta2[1]*x1 + beta2[2]*x2 + beta2[3]*B + beta2[4]*C + 
                                     beta2[5]*D1 + beta2[6]*D2 + beta2[7]*D3, lin_pred)) %>%
  dplyr::mutate(lin_pred = if_else(latent == 3,
                                   group_intercepts_3[level] +  beta3[1]*x1 + beta3[2]*x2 + beta3[3]*B + beta3[4]*C + 
                                     beta3[5]*D1 + beta3[6]*D2 + beta3[7]*D3, lin_pred ))


y <- vector("numeric",dim(data)[1])
for(i in 1:dim(data)[1]){
  y[i] <- rbinom(1, size = 1, prob = plogis(data$lin_pred[i]))
}
data <- cbind(y,data)
data <- data[, !names(data) %in% 'lin_pred']
data = onehot_to_cat(data, list(c("A", "B", "C")), "V1")

test_data <- data[sample(1:nrow(data), 200, replace = FALSE), ]

saveRDS(test_data, file=filename_data)
  
