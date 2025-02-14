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

  # Simulate data
  data <- as.data.frame(rbind(mvrnorm(n = w[1]*n_obs, mu = mean.1, Sigma = cov.1), mvrnorm(n = w[2]*n_obs, mu = mean.2, Sigma = cov.2) ))
  data <- cbind(rep(1:n_groups, each = n_per_group)[sample(n_obs)],
                c(rep(1, length.out = w[1]*n_obs),rep(2, length.out = w[2]*n_obs)), data)
  colnames(data) <- c('level','latent','x1','x2')
  
  result1 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = cat_probs1),
                   rmultinom(w[2]*n_obs, size = 1, prob = cat_probs2))
  
  V <- data.frame( v1 = data.frame(t(result1)) )
  colnames(V) <- c("A","B","C")
  data <- cbind(data,V)
  
  dich.1 <- as.data.frame(IsingSampler(w[1]*n_obs,int1,thres1,beta = 1))
  dich.2 <- as.data.frame(IsingSampler(w[2]*n_obs,int2,thres2, beta = 1))
  dich <- rbind(dich.1,dich.2)
  colnames(dich) <- c("D1","D2","D3")
  data <- cbind(data,dich)

  data <- data %>%
    dplyr::mutate(lin_pred = 0) %>%
    dplyr::mutate(lin_pred = if_else(latent == 1,
                                     group_intercepts_1[level] + beta1[1]*x1 + beta1[2]*x2 + beta1[3]*B + beta1[4]*C + 
                                       beta1[5]*D1 + beta1[6]*D2 + beta1[7]*D3,lin_pred)) %>%
    dplyr::mutate(lin_pred = if_else(latent == 2, 
                                     group_intercepts_2[level] + beta2[1]*x1 + beta2[2]*x2 + beta2[3]*B + beta2[4]*C + 
                                       beta2[5]*D1 + beta2[6]*D2 + beta2[7]*D3, lin_pred)) 

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
  
  saveRDS(data, file=filename_data)

}


# Get test data
set.seed(3792)
filename_data <- paste("sim_data/testdata", ".rds", sep = "")

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

# Simulate data
data <- as.data.frame(rbind(mvrnorm(n = w[1]*n_obs, mu = mean.1, Sigma = cov.1), mvrnorm(n = w[2]*n_obs, mu = mean.2, Sigma = cov.2) ))
data <- cbind(rep(1:n_groups, each = n_per_group)[sample(n_obs)],
              c(rep(1, length.out = w[1]*n_obs),rep(2, length.out = w[2]*n_obs)), data)
colnames(data) <- c('level','latent','x1','x2')

result1 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = cat_probs1),
                 rmultinom(w[2]*n_obs, size = 1, prob = cat_probs2))

V <- data.frame( v1 = data.frame(t(result1)) )
colnames(V) <- c("A","B","C")
data <- cbind(data,V)

dich.1 <- as.data.frame(IsingSampler(w[1]*n_obs,int1,thres1,beta = 1))
dich.2 <- as.data.frame(IsingSampler(w[2]*n_obs,int2,thres2, beta = 1))
dich <- rbind(dich.1,dich.2)
colnames(dich) <- c("D1","D2","D3")
data <- cbind(data,dich)

data <- data %>%
  dplyr::mutate(lin_pred = 0) %>%
  dplyr::mutate(lin_pred = if_else(latent == 1,
                                   group_intercepts_1[level] + beta1[1]*x1 + beta1[2]*x2 + beta1[3]*B + beta1[4]*C + 
                                     beta1[5]*D1 + beta1[6]*D2 + beta1[7]*D3,lin_pred)) %>%
  dplyr::mutate(lin_pred = if_else(latent == 2, 
                                   group_intercepts_2[level] + beta2[1]*x1 + beta2[2]*x2 + beta2[3]*B + beta2[4]*C + 
                                     beta2[5]*D1 + beta2[6]*D2 + beta2[7]*D3, lin_pred)) 

y <- vector("numeric",dim(data)[1])
for(i in 1:dim(data)[1]){
  y[i] <- rbinom(1, size = 1, prob = plogis(data$lin_pred[i]))
}
data <- cbind(y,data)
data <- data[, !names(data) %in% 'lin_pred']
data = onehot_to_cat(data, list(c("A", "B", "C")), "V1")

test_data <- data[sample(1:nrow(data), 200, replace = FALSE), ]

saveRDS(test_data, file=filename_data)

