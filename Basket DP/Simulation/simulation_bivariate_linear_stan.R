# Bivariate linear in Stan

library(mvtnorm)
library(systemfit)
library(rstan)
set.seed(666)
plot.wd	<- "/kellogg/users/marketing/2661703/Exercise"

# Simulate data 
beta<- cbind(c(.8, -1, .5), c(.5, -.5, 1))
rho	<- .5
N	<- 1000				# Number of observations
D	<- 2				# Number of dimensions of responses
K	<- 3				# Number of covariates
Sigma<- matrix(c(1, rho, rho, 1), D, D)
X	<- cbind(1, matrix(runif(N*(K-1), -1, 1), N, K-1))
mu	<- X %*% beta
ystar	<- t(sapply(1:N, function(i) rmvnorm(1, mu[i,], Sigma)))
y	<- 1*(ystar>0)
colMeans(y)

# SUR model 
mydata	<- data.frame(y= y, X[,-1])
f1	<- systemfit(list(eq1 = y.1 ~ X1 + X2, eq2 = y.2 ~ X1 + X2), data = mydata, method = "SUR")
summary(f1)

# Stan source code
src_code <- "
data{
	int<lower = 1>	N; 
	int<lower = 1>	D; 
	int<lower = 1>	K; 
	vector[K]		X[N];
	vector[D] 		y[N];
	cov_matrix[D]	ID;
}
parameters{
	vector[K]	beta[D];
	cov_matrix[D] Sigma;
}
model{
	matrix[N,D] Mu; 
	
	# Prior distribution
	for(i in 1:D){
		beta[i] ~ normal(0, 10);
	}
	Sigma ~ inv_wishart(D + 1, ID);
	# Linear response
	for(i in 1:N){
		for(j in 1:D){
			Mu[i,j]	<- dot_product(X[i], beta[j]); 
		}
	}
	
	# The distribution of discrete outcome
	for(i in 1:N){
			y[i] ~ multi_normal(row(Mu, i)', Sigma); 
	}
}
"
mydata		<- list(N = N, D = D, K = K, X = X, y = y, ID = diag(D))
save_par	<- c("beta", "Sigma")
my_init		<- function(){
	list(beta = t(beta), Sigma = Sigma, z = ystar + matrix(rnorm(N*D), N, D) )
}
m	<- stan_model(model_code = src_code)
f2 	<- optimizing(m, data = mydata, hessian = TRUE)

# See if SUR and Bivariate normal generates the same results 
par	<- f2$par
sel	<- grep("beta", names(par))
par[sel]
coef(f1)

# Covariance matrix
sel	<- grep("Sigma", names(par))
par[sel]
f1$residCov
