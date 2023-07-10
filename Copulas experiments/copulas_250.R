rm(list = ls())
library(MASS)
library(compositions)
library(EasyABC)
library(copula)

set.seed(1)

df <-read.table(file.choose(), header=FALSE, strip.white=TRUE)
attach(df)
df <- cbind(V1,V2)

# Data
n = 250
theta = 0.3
cop.1 <- mvdc(gumbelCopula(2), margins = c("norm", "norm"), paramMargins = list(list(0, 1), list(0, 1)))
cop.2 <- mvdc(claytonCopula(2.13), margins = c("norm", "norm"), paramMargins = list(list(0, 1), list(0, 1)))
tmp <- runif(n, 0, 1)
X <- matrix(0, n, 2)
for (i in 1:n) if (tmp[i] < theta) X[i, ] <- rMvdc(1, cop.1) else X[i, ] <- rMvdc(1, cop.2)

r_mat <- c("true posterior prob", 
			"t1", "t2", "t3", "t4", "t5", "et",
			"T1", "T2", "T3", "T4", "T5", "eT",
			"M1", "M2", "M3", "M4", "M5", "eM", 
			"C1", "C2", "C3", "C4", "C5", "eC")

for (I in 1:dim(df)[1]){
cat("\n", df[I, ], "\n")

# The intervals
int <- df[I, ]

# Define the Likelihood
f <- function(t, x) {
	n <- dim(x)[1]
	result <- 1
	for (i in 1:n) result <- result * (t*dMvdc(x[i, ], cop.1) + (1 - t)*dMvdc(x[i, ], cop.2))
	result * dunif(t, 0, 1)
}
g <- function(t) f(t, X)

# True posterior probability
output <- integrate(g, int[1], int[2])$value
denominator <- integrate(g, 1e-4, 1 - 1e-4)$value

r_vec <- (output/denominator)
print(r_vec)

# directional vectors
n_u <- 10
u <- matrix(0, n_u, 2)
for (i in 1:n_u) {
	u[i, 1] <- cos(i*2*pi/n_u)
	u[i, 2] <- sin(i*2*pi/n_u)
}

#1. TDC 0.05
alpha <- c(.05)
n_alpha <- length(alpha)

toy_model <- function(t){
	tmp <- runif(n, 0, 1)
	x <- matrix(0, n, 2)
	for (i in 1:n) if (tmp[i] < t) x[i, ] <- rMvdc(1, cop.1) else x[i, ] <- rMvdc(1, cop.2)
	Q <- replicate(n_u*n_alpha, NA)
	for (k in 1:n_alpha) 
		for (i in 1:n_u) Q[i+(k-1)*n_u] <- quantile(x%*%u[i, ], alpha[k])
	return(Q)
}

# Compute Summary Statistics of the observations
Q <- replicate(n_u*n_alpha, NA)
for (k in 1:n_alpha) 
	for (i in 1:n_u) Q[i+(k-1)*n_u] <- quantile(X%*%u[i, ], alpha[k])
sum_stat_obs=Q

# Approximate Bayesian Computation
toy_prior=list(c("unif",0,1))

N = 10^4
p = 0.01
K = N*p

# Calculate the portion
n_mean <- 5
R <- rep(NA, n_mean)
for (i in 1:n_mean) {
	ABC_rej<-ABC_rejection(model=toy_model, prior=toy_prior, nb_simul=N,
					summary_stat_target=sum_stat_obs, tol=p)
	TF <- rep(0, K)
	for (j in 1:K) {
		tmp <- ((ABC_rej$param[j] > int[1]) & (ABC_rej$param[j] < int[2]))
		TF[j] <- all(tmp == 1)
	}
	R[i] <- sum(TF)/(N*p)
}
r_vec <- c(r_vec, R)
r_vec <- c(r_vec, mean(abs(R-output/denominator)))
print(r_vec)

#2. TDC 0.4
alpha <- c(.4)
n_alpha <- length(alpha)

toy_model <- function(t){
	tmp <- runif(n, 0, 1)
	x <- matrix(0, n, 2)
	for (i in 1:n) if (tmp[i] < t) x[i, ] <- rMvdc(1, cop.1) else x[i, ] <- rMvdc(1, cop.2)
	Q <- replicate(n_u*n_alpha, NA)
	for (k in 1:n_alpha) 
		for (i in 1:n_u) Q[i+(k-1)*n_u] <- quantile(x%*%u[i, ], alpha[k])
	return(Q)
}

# Compute Summary Statistics of the observations
Q <- replicate(n_u*n_alpha, NA)
for (k in 1:n_alpha) 
	for (i in 1:n_u) Q[i+(k-1)*n_u] <- quantile(X%*%u[i, ], alpha[k])
sum_stat_obs=Q

# Approximate Bayesian Computation
toy_prior=list(c("unif",0,1))

N = 10^4
p = 0.01
K = N*p

# Calculate the portion
n_mean <- 5
R <- rep(NA, n_mean)
for (i in 1:n_mean) {
	ABC_rej<-ABC_rejection(model=toy_model, prior=toy_prior, nb_simul=N,
					summary_stat_target=sum_stat_obs, tol=p)
	TF <- rep(0, K)
	for (j in 1:K) {
		tmp <- ((ABC_rej$param[j] > int[1]) & (ABC_rej$param[j] < int[2]))
		TF[j] <- all(tmp == 1)
	}
	R[i] <- sum(TF)/(N*p)
}
r_vec <- c(r_vec, R)
r_vec <- c(r_vec, mean(abs(R-output/denominator)))
print(r_vec)

#3. Multiple TDCs
alpha <- c(.05, .1, .2, .4)
n_alpha <- length(alpha)

toy_model <- function(t){
	tmp <- runif(n, 0, 1)
	x <- matrix(0, n, 2)
	for (i in 1:n) if (tmp[i] < t) x[i, ] <- rMvdc(1, cop.1) else x[i, ] <- rMvdc(1, cop.2)
	Q <- replicate(n_u*n_alpha, NA)
	for (k in 1:n_alpha) 
		for (i in 1:n_u) Q[i+(k-1)*n_u] <- quantile(x%*%u[i, ], alpha[k])
	return(Q)
}

# Compute Summary Statistics of the observations
Q <- replicate(n_u*n_alpha, NA)
for (k in 1:n_alpha) 
	for (i in 1:n_u) Q[i+(k-1)*n_u] <- quantile(X%*%u[i, ], alpha[k])
sum_stat_obs=Q

# Approximate Bayesian Computation
toy_prior=list(c("unif",0,1))

N = 10^4
p = 0.01
K = N*p

# Calculate the portion
n_mean <- 5
R <- rep(NA, n_mean)
for (i in 1:n_mean) {
	ABC_rej<-ABC_rejection(model=toy_model, prior=toy_prior, nb_simul=N,
					summary_stat_target=sum_stat_obs, tol=p)
	TF <- rep(0, K)
	for (j in 1:K) {
		tmp <- ((ABC_rej$param[j] > int[1]) & (ABC_rej$param[j] < int[2]))
		TF[j] <- all(tmp == 1)
	}
	R[i] <- sum(TF)/(N*p)
}
r_vec <- c(r_vec, R)
r_vec <- c(r_vec, mean(abs(R-output/denominator)))
print(r_vec)

#4. Common Statistics
toy_model <- function(t){
	tmp <- runif(n, 0, 1)
	x <- matrix(0, n, 2)
	for (i in 1:n) if (tmp[i] < t) x[i, ] <- rMvdc(1, cop.1) else x[i, ] <- rMvdc(1, cop.2)
	C <- var(x)
	S <- c(colMeans(x), C[1,1], C[1,2], C[2,2])
	return(S)
}

# Compute Summary Statistics of the observations
C <- var(X)
S <- c(colMeans(X), C[1,1], C[1,2], C[2,2])
sum_stat_obs=S

# Approximate Bayesian Computation
toy_prior=list(c("unif",0,1))

N = 10^4
p = 0.01
K = N*p

# Calculate the portion
n_mean <- 5
R <- rep(NA, n_mean)
for (i in 1:n_mean) {
	ABC_rej<-ABC_rejection(model=toy_model, prior=toy_prior, nb_simul=N,
					summary_stat_target=sum_stat_obs, tol=p)
	TF <- rep(0, K)
	for (j in 1:K) {
		tmp <- ((ABC_rej$param[j] > int[1]) & (ABC_rej$param[j] < int[2]))
		TF[j] <- all(tmp == 1)
	}
	R[i] <- sum(TF)/(N*p)
}
r_vec <- c(r_vec, R)
r_vec <- c(r_vec, mean(abs(R-output/denominator)))
print(r_vec)

r_mat <- rbind(r_mat, r_vec)
}
write.csv(r_mat, "copulas_250_1.csv")