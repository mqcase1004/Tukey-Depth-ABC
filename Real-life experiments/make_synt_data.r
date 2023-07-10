rm(list = ls())

F = function(time, t, ts) {
	val = 0;
	if(time < t) val = 0.0
	else if( (t<time) & (ts>time) ) val = time - t
	else val = ts - t
	return(val)
}

D = 40 # Number of days with observations
N = 20 # Number of events per day
T = 144 # Number of observations per day (288)
time = seq(0, 1, length.out = T+1)

shape1 = 6
shape2 = 4
par(mfrow = c(3,2))
plot(time, dbeta(time, shape1, shape2), type = "l", col = "red",
xlab = "arrival time", ylab = "density", ylim = c(0,5.5),
xaxt = "n",
main = expression(paste("Beta distribution, ", phi, " = 6, ", beta, " = 4")))
axis(1, at = seq(2,22,4)/24, labels = c("2 am", "6 am", "10 am", "2 pm", "6 pm", "10 pm"))

# Generate synthetic data
t = matrix(nrow = D, ncol = N)
for(d in 1:D) {
	t[d,] = rbeta(n = N, shape1 = shape1, shape2 = shape2)
}

s = matrix(nrow = D, ncol = N)
rateexp = 10
for(d in 1:D) {
	s[d,] = rexp(n = N, rate = rateexp)
}

plot(seq(0, 10/rateexp, 1e-2), dexp(seq(0, 10/rateexp, 1e-2), rate = rateexp), type = "l",
xlab = "processing time", ylab = "density", col = "red",
main = expression(paste("Exponential distribution, ", lambda, " = 10")))


y = matrix(nrow = D, ncol = T)
for(d in 1:D) {
	for(j in 1:T) {
		m = 0.0
		for(i in 1:N) {
			m = m + (F(time[j+1], t[d,i], t[d,i]+s[d,i]) - F(time[j], t[d,i], t[d,i]+s[d,i]))/(time[j+1] - time[j])
		}
	y[d,j] = m # no additinal random noise
	}
}

plot(c(t[1,1], t[1,1]+s[1,1]), c(1, 1), xlim = c(0,1), ylim = c(0,N), type = "l", lwd = 4)
for(i in 2:N) {
	points(c(t[1,i], t[1,i]+s[1,i]), i*c(1, 1), xlim = c(0,1), ylim = c(0,N), type = "l", lwd = 4)
}

plot(c(t[2,1], t[2,1]+s[2,1]), c(1, 1), xlim = c(0,1), ylim = c(0,N), type = "l", lwd = 4)
for(i in 2:N) {
	points(c(t[2,i], t[2,i]+s[2,i]), i*c(1, 1), xlim = c(0,1), ylim = c(0,N), type = "l", lwd = 4)
}

plot(y[1,], type = "l", col = "blue")
plot(y[2,], type = "l", col = "blue")

S = matrix(nrow = 3, ncol = D)
for (d in 1:D) {
	S[1, d] = mean(y[d, ])
	S[2, d] = sum((1:T)*(y[d, ]/sum(y[d, ])))
	S[3, d] = sum((1:T)^2*(y[d, ]/sum(y[d, ]))) - (S[2, d])^2
}

# Stats <- c(S[1, ], S[2, ], S[3, ])