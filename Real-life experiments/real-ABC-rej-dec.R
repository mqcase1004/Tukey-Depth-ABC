rm(list = ls())
library(MASS)
library(EasyABC)

F = function(time, t, ts) {
	val = 0;
	if(time < t) val = 0
	else if( (t<time) & (ts>time) ) val = time - t
	else val = ts - t
	return(val)
}

set.seed(1)
D = 40 # Number of days with observations
N = 20 # Number of events per day
T = 96 # Number of observations per day (288)
time = seq(0, 1, length.out = T+1)

# directional vectors
n_u <- 10
u <- mvrnorm(n_u, rep(0, 3), diag(3))
for (i in 1:n_u) u[i, ] <- u[i, ]/norm(u[i, ], type = "2")

# Tukey Depth Contour 0.05
alpha <- .05

# Generate synthetic data
model <- function(w){
	t <- matrix(nrow = D, ncol = N)
	for(d in 1:D) t[d,] <- rbeta(n = N, shape1 = w[1], shape2 = w[2])
	s <- matrix(nrow = D, ncol = N)
	for(d in 1:D) s[d,] <- rexp(n = N, rate = w[3])
	y <- matrix(nrow = D, ncol = T)
	for(d in 1:D) 
		for(j in 1:T) {
			m <- 0
			for(i in 1:N) 
				m <- m + (F(time[j+1], t[d,i], t[d,i]+s[d,i]) - F(time[j], t[d,i], t[d,i]+s[d,i]))/(time[j+1] - time[j])
			y[d,j] <- m # no additinal random noise
		}
	Stats <- matrix(nrow = D, ncol = 3)
	Stats[, 1] <- rowMeans(y)
	tmp <- y/rowSums(y)
	Stats[, 2] <- tmp%*%(1:T)
	Stats[, 3] <- tmp%*%((1:T)^2) - (Stats[, 2])^2
	rescale = c(1, 60, 200)
	Stats <- t(t(Stats)/rescale)
	Q <- apply(Stats%*%t(u), 2, quantile, probs = alpha)
	return(Q)
}

# Compute Summary Statistics of the observations
Q <- model(c(6, 4, 10))
sum_stat_obs=Q

# Approximate Bayesian Computation
prior=list(c("unif",3,8), c("unif",2,7), c("unif",3,40))

N_abc = 1e4
p = 1e-2
K = N_abc*p

start_time <- Sys.time()
ABC_rej<-ABC_rejection(model, prior, nb_simul=N_abc, summary_stat_target=sum_stat_obs, tol=p)
end_time <- Sys.time()
end_time - start_time

result <- ABC_rej$param
result
write.csv(result, "real-result-dec10-40-20-96.csv", row.names=FALSE)
