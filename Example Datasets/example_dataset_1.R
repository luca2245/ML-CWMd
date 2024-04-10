set.seed(444)

n_obs <- 2000

n_groups <- 10 # number of levels(hospitals)
n_per_group <- 200 # number of observations per level(hospital)
group_ids <- rep(1:n_groups, each = n_per_group)

w <- c(0.20,0.30,0.50) # proportion of observations per cluster

# Random Intercept
group_intercepts_1 <- rnorm(n_groups, mean = 0, sd = 2)
group_intercepts_2 <- rnorm(n_groups, mean = 0, sd = 2)
group_intercepts_3 <- rnorm(n_groups, mean = 0, sd = 2)

# MU and SIGMA Parameters
cov.1 <- matrix(c(0.7,0.5,0.5,3),nrow = 2)
mean.1 <- c(2.04798648, 0.12524515)

cov.2 <- matrix(c(2,-1,-1,3),nrow = 2)
mean.2 <- c(5.057435, 4.839869)

cov.3 <- matrix(c(3,1,1,2),nrow = 2)
mean.3 <- c(4.219037, -4.509953)

# LAMBDA Parameters
probs1 <- c(0.74, 0.18, 0.08)

s_probs <- c(0.08, 0.74, 0.18)

probs1.1 <- c(0.30, 0.51, 0.19)

probs2 <- c(0.51, 0.49)

probs2.1 <- c(0.47, 0.53)

#INTERACTION AND THRESHOLD Parameters
N <- 3

int1 <- matrix(c(0, 0.21, -1.10, 0.21, 0, 0, -1.10, 0, 0), nrow = 3, ncol = 3)
thres1 <- c(0.1093606, 0.3826030, -0.4929232)


int2 <- matrix(c(0, -0.73, 0.83, -0.73,  0, 0, 0.83,  0, 0), nrow = 3, ncol = 3)
thres2 <- c(-0.1450787,  0.8800758, -0.1818312)

int3 <- matrix(c(0, -4.15, 2.11, -4.15,  
                 0,  1.14, 2.11,  1.14,  0), nrow = 3, ncol = 3)
thres3 <- c(0.73583084, -0.23133994,  0.01153844)

# Fixed effects
beta1 <- c(-0.52089065,  0.08335639, -1.30575016,  0.22371244,  5.33094751,  2.75188423, -2.29292748,  0.93377294) 

beta2 <- c(-0.06715653,  0.79119887, -0.45963063,  0.25307799, -3.88824880, -0.63264122,  0.62944698, -1.51371392)

beta3 <- c(-0.4216015, -0.3076446, -1.3255474, -0.6001705, -4.1797789,  4.8872690,  3.3393496, -0.4575995)

data <- as.data.frame(rbind(mvrnorm(n = w[1]*n_obs, mu = mean.1, Sigma = cov.1), mvrnorm(n = w[2]*n_obs, mu = mean.2, Sigma = cov.2), mvrnorm(n = w[3]*n_obs, mu = mean.3, Sigma = cov.3)))
data <- cbind(rep(1:n_groups, each = n_per_group)[sample(n_obs)],
              c(rep(1, length.out = w[1]*n_obs),rep(2, length.out = w[2]*n_obs),
                rep(3, length.out = w[3]*n_obs)), data)
colnames(data) <- c('level','latent','x1','x2')

result1 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = probs1),
                 rmultinom(w[2]*n_obs, size = 1, prob = s_probs),
                 rmultinom(w[3]*n_obs, size = 1, prob = probs1.1))

result2 <- cbind(rmultinom(w[1]*n_obs, size = 1, prob = probs2),
                 rmultinom(w[2]*n_obs, size = 1, prob = 1 - probs2),
                 rmultinom(w[3]*n_obs, size = 1, prob = probs2.1))

V <- data.frame(v1 = data.frame(t(result1)), v2 = data.frame(t(result2)))
colnames(V) <- c("A","B","C","a","b")
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
                                   group_intercepts_1[level] + beta1[1]*x1 + beta1[2]*x2 + beta1[3]*A + beta1[4]*B + beta1[5]*a + 
                                     beta1[6]*D1 + beta1[7]*D2 + beta1[8]*D3,lin_pred)) %>%
  dplyr::mutate(lin_pred = if_else(latent == 2, group_intercepts_2[level] + beta2[1]*x1 + beta2[2]*x2 + beta2[3]*A +              beta2[4]*B + beta2[5]*a + 
                                     beta2[6]*D1 + beta2[7]*D2 + beta2[8]*D3, lin_pred)) %>%
  dplyr::mutate(lin_pred = if_else(latent == 3,
                                   group_intercepts_3[level] +  beta3[1]*x1 + beta3[2]*x2 + beta3[3]*A + beta3[4]*B + beta3[5]*a + 
                                     beta3[6]*D1 + beta3[7]*D2 + beta3[8]*D3, lin_pred ))


y <- vector("numeric",dim(data)[1])
for(i in 1:dim(data)[1]){
  y[i] <- rbinom(1, size = 1, prob = plogis(data$lin_pred[i]))
}
data <- cbind(y,data)

data <- data[sample(nrow(data)),]

check.1 <- sum(data$y[which(data$latent == 1)])/(w[1]*n_obs)
check.2 <- sum(data$y[which(data$latent == 2)])/(w[2]*n_obs)
check.3 <- sum(data$y[which(data$latent == 3)])/(w[3]*n_obs)


U <- data[,4:5]
V1 <- data[,6:8]
V2 <- data[,9:10]
V <- list(V1,V2)
D <- data[,11:13]
Y <- data[,1]

formula <- y ~ x1 + x2 +  A + B + a + D1 + D2 + D3 + (1|level) - 1

v <- 2 
d <- 3 
