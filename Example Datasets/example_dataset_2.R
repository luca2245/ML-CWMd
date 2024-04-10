# Create Example Dataset
n_obs <- 5000
# Set levels(hospitals)
n_groups <- 10 # number of levels
n_per_group <- 500 # number of observations per level
group_ids <- rep(1:n_groups, each = n_per_group)

w <- c(0.20,0.30,0.50)

#random intercept
set.seed(2023)
group_intercepts_1 <- rnorm(n_groups, mean = -1, sd = 2)
group_intercepts_2 <- rnorm(n_groups, mean = 0, sd = 2)
group_intercepts_3 <- rnorm(n_groups, mean = 1, sd = 2)

mean.1 <- c(0,0)
cov.1 <- matrix(c(0.7,0.5,0.5,3),nrow = 2)

mean.2 <- c(6,6)
cov.2 <- matrix(c(2,-1,-1,3),nrow = 2)

mean.3 <- c(3.5,-4)
cov.3 <- matrix(c(3,1,1,2),nrow = 2)

data <- as.data.frame(rbind(mvrnorm(n = w[1]*n_obs, mu = mean.1, Sigma = cov.1), mvrnorm(n = w[2]*n_obs, mu = mean.2, Sigma = cov.2), mvrnorm(n = w[3]*n_obs, mu = mean.3, Sigma = cov.3)))
data <- cbind(rep(1:n_groups, each = n_per_group)[sample(n_obs)],
              c(rep(1, length.out = w[1]*n_obs),rep(2, length.out = w[2]*n_obs),
                rep(3, length.out = w[3]*n_obs)), data)
colnames(data) <- c('level','latent','x1','x2')

# Categorical:
n_obs <- 5000

probs1 <- c(0.2, 0.1, 0.7)
probs2 <- c(0.4, 0.6)
probs3 <- c(0.2,0.8)


result1 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = probs1),
                 rmultinom(w[2]*n_obs, size = 1, prob = c(0.35,0.45,0.2)),
                 rmultinom(w[3]*n_obs, size = 1, prob = c(0.4,0.3,0.3)))

result2 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = probs2),
                 rmultinom(w[2]*n_obs, size = 1, prob = 1 - probs2),
                 rmultinom(w[3]*n_obs, size = 1, prob = c(0.5,0.5)))

result3 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = probs3),
                 rmultinom(w[2]*n_obs, size = 1, prob = 1 - probs3),
                 rmultinom(w[3]*n_obs, size = 1, prob = c(0.3,0.7)))

V <- data.frame(v1 = data.frame(t(result1)), v2 = data.frame(t(result2)), v3 = data.frame(t(result3)))
colnames(V) <- c("A","B","C","a","b","D","E")
data <- cbind(data,V)

# Generating random dependent dichotomic variables using IsingSampler:

set.seed(203)
N <- 5
int1 <- matrix(sample(0:1,N^2,TRUE,prob = c(0.4, 0.6)),N,N) * rnorm(N^2)
int1 <- int1 + t(int1)
diag(int1) <- 0
thres1 <- rnorm(5,0,0.5)



set.seed(145)
int2 <- matrix(sample(0:1,N^2,TRUE,prob = c(0.7, 0.3)),N,N) * rnorm(N^2,0,1.5)
int2 <- int2 + t(int2)
diag(int2) <- 0
thres2 <- rnorm(5,0,0.5)



set.seed(180)
int3 <- matrix(sample(0:1,N^2,TRUE,prob = c(0.6, 0.4)),N,N) * rnorm(N^2,0,2)
int3 <- int3 + t(int3)
diag(int3) <- 0
thres3 <- rnorm(5,0,0.5)



dich.1 <- as.data.frame(IsingSampler(w[1]*n_obs,int1,thres1,beta = 1))
dich.2 <- as.data.frame(IsingSampler(w[2]*n_obs,int2,thres2, beta = 1))
dich.3 <- as.data.frame(IsingSampler(w[3]*n_obs,int3,thres3, beta = 1))
dich <- rbind(dich.1,dich.2,dich.3)
colnames(dich) <- c("D1","D2","D3","D4","D5")
data <- cbind(data,dich)


#Linear Predictor:

data <- data %>%
  dplyr::mutate(lin_pred = 0) %>%
  dplyr::mutate(lin_pred = if_else(latent == 1,
                                   group_intercepts_1[level] + 2*x1 + -0.4*x2 + 1.5*A + 0.7*B + -0.8*a + -1.1*D + 
                                     0.3*D1 -0.3*D2 -2*D3 + 3.6*D4 -0.2*D5,lin_pred)) %>%
  dplyr::mutate(lin_pred = if_else(latent == 2, group_intercepts_2[level] + 0.8*x1 - 0.6*x2 + -0.9*A + 
                                     1.5*B + -0.5*a - 1.4*D -0.3*D1 -0.3*D2 +3*D3 -1.7*D4 -0.05*D5, lin_pred)) %>%
  dplyr::mutate(lin_pred = if_else(latent == 3,
                                   group_intercepts_3[level] + -0.8*x1 + 0.4*x2 + -0.9*A + 1.5*B + 1.2*a + +1.4*D +
                                     0.7*D1 -0.8*D2 +1.88*D3 -0.5*D4 + 1.5*D5, lin_pred ))




y <- vector("numeric",dim(data)[1])
for(i in 1:dim(data)[1]){
  y[i] <- rbinom(1, size = 1, prob = plogis(data$lin_pred[i]))
}
data <- cbind(y,data)


data <- data[sample(nrow(data)),]

U <- data[,4:5]
V1 <- data[,6:8]
V2 <- data[,9:10]
V3 <- data[,11:12]
V <- list(V1,V2,V3)
D <- data[,13:17]
Y <- data[,1]

formula <- y ~ x1 + x2 +  A + B + a + D + D1 + D2 + D3 + D4 + D5 + (1|level)

v <- 3 # number of categorical covariates
d <- 5 #number of dependent dichotomical variables
