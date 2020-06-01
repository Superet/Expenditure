# Bivariate probit in Stan

library(mvtnorm)
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

# Stan source code
src_code <- "
data{
	int<lower = 1>	N; 
	int<lower = 1>	D; 
	int<lower = 1>	K; 
	vector[K]		X[N];
	int 			y[N,D];
	cov_matrix[D]	ID;
}
parameters{
	vector[K]	beta[D];
	corr_matrix[D] Sigma;
	vector[D]	z[N];
}
transformed parameters{
	real        scale; 
	vector[K]	beta_adj[D];
	corr_matrix[D]	Sigma_adj;
	
	scale 	<- sqrt(2); 
	for(i in 1:D){
		beta_adj[i]	<- beta[i]/scale; 
	}
	Sigma_adj	<- Sigma; 
	Sigma_adj[1,2]	<- Sigma[1,2]/scale; 
	Sigma_adj[2,1]	<- Sigma[2,1]/scale;
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
		z[i] ~ multi_normal(row(Mu, i)', Sigma); 
	}
	for(i in 1:N){
		for(j in 1:D){
			y[i,j] ~ bernoulli(Phi(z[i,j])); 
		}
	}
}
"

mydata		<- list(N = N, D = D, K = K, X = X, y = y, ID = diag(D))
save_par	<- c("beta", "Sigma", "beta_adj","Sigma_adj")
my_init		<- function(){
	list(beta = t(beta), Sigma = Sigma, z = ystar + matrix(rnorm(N*D), N, D) )
}
prc			<- proc.time()
fit 		<- stan(model_code = src_code, data = mydata, pars = save_par, init = my_init,
					chains = 4, iter = 120, warmup = 20, thin = 4)
use.time	<- proc.time() - prc
cat("Rstan took", use.time[3]/60, "min.\n")					
print(fit)

pdf(paste(plot.wd, "/graph_stanfit.pdf",sep=""), width = 10, height = 8)
plot(fit)
dev.off() 

# MCMCglmm
library(MCMCglmm)
mydat 	<- data.frame(y = y, X = X[,-1])
m1		<- MCMCglmm(cbind(y.1, y.2) ~ trait:X.1 + trait:X.2 + trait -1, rcov = ~ us(trait):units, 
					family = c("categorical","categorical"), data = mydat,  nitt = 11000, burnin = 1000, thin= 10)
colMeans(m1$Sol)
colMeans(m1$VCV)
