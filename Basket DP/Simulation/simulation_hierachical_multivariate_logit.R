# Simulation of hierachical multivariate logit model 

# Latent l(N*D)
# l = X*beta + beta0 + epsilon
# 


library(mvtnorm)
library(rstan)
library(doParallel)
library(parallel)
library(foreach)
library(MCMCglmm)
set.seed(666)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
plot.wd	<- "/kellogg/users/marketing/2661703/Exercise"

options(error = quote({dump.frames(to.file = TRUE)}))

# Simulate data 
beta	<- cbind(c(.8, -1, .5), c(.5, -.5, 1), c(0, .2, .3))
rho12	<- .3
rho13	<- -.5
rho23	<- .4
N		<- 1000				# Number of observations
D		<- 3				# Number of dimensions of responses
K		<- 3				# Number of covariates
J		<- 20				# Number of groups
omega	<- c(1, .6, 1.2)
beta0	<- sapply(1:D, function(i) rnorm(J, 0, omega[i]))
Sigma	<- diag(D)
Sigma[lower.tri(Sigma)]	<- c(rho12, rho13, rho23)
Sigma[upper.tri(Sigma)]	<- c(rho12, rho13, rho23)
X		<- cbind(1, matrix(runif(N*(K-1), -1, 1), N, K-1))
id.idx	<- sample(1:J, N, replace = T)
mu		<- X %*% beta + beta0[id.idx,]
ystar	<- t(sapply(1:N, function(i) rmvnorm(1, mu[i,], Sigma)))
p		<- exp(ystar)/(1 + exp(ystar))
y		<- apply(p, 2, function(x) rbinom(N, 1, x))
colMeans(y)

# Stan source code
src_code <- "
data{
	int<lower = 1>	N; 
	int<lower = 1>	D; 
	int<lower = 1>	K; 
	int<lower = 1>	J; 
	vector[K]		X[N];
	int 			y[N,D];
	cov_matrix[D]	ID;
	int<lower = 1> 	idx[N]; 
}
parameters{
	vector[K]	beta[D];
	vector[J]	beta0[D];
	vector<lower = 0>[D]	omega;
	corr_matrix[D] Sigma;
	vector[D]	z[N];
}
model{
	matrix[N,D] Mu; 
	
	# Prior distribution
	for(i in 1:D){
		beta[i] ~ normal(0, 10);
	}
	omega ~ uniform(0, 10); 
	for(i in 1:D){
		beta0[i] ~ normal(0, omega[i]); 
	}
	Sigma ~ inv_wishart(D + 1, ID);
	
	# Linear response
	for(i in 1:N){
		for(j in 1:D){
			Mu[i,j]	<- dot_product(X[i], beta[j]) + beta0[j,idx[i]] ; 
		}
	}
	
	# The distribution of discrete outcome
	for(i in 1:N){
		z[i] ~ multi_normal(row(Mu, i)', Sigma); 
	}
	for(i in 1:N){
		for(j in 1:D){
			y[i,j] ~ bernoulli_logit(z[i,j]); 
		}
	}
}
"

mydata		<- list(N = N, D = D, K = K, J = J, X = X, y = y, ID = diag(D), idx = id.idx, zz = rep(0, D))
save_par	<- c("beta", "Sigma","omega", "beta0")
init_ls		<- list(beta = t(as.matrix(beta)), beta0 = t(as.matrix(beta0)), omega = omega, Sigma = Sigma, z = ystar + matrix(rnorm(N*D), N, D) )
my_init		<- function(){return(init_ls)}
foo			<- stan(model_code = src_code, data = mydata, pars = save_par, init = my_init, chains = 1, iter = 1)

CL = makeCluster(3)
clusterExport(cl = CL, c("foo", "src_code", "mydata", "save_par", "my_init", "init_ls")) 
prt		<- proc.time()
fitls <- parLapply(CL, 1:3, fun = function(cid) {  
  	require(rstan)
	stan(fit = foo, data = mydata, pars = save_par, init = my_init,
					chains = 1, iter = 120, warmup = 20, thin = 4, chain_id = cid)
})
use.time1	<- proc.time() - prt
stopCluster(CL)
cat("Stan sampling takes", use.time1[3]/60,"min.\n")

fit		<- sflist2stanfit(fitls)
print(fit)

# MCMCglmm #
mydat 	<- data.frame(y = y, X = X[,-1], id = id.idx)
names(mydat)	<- gsub("\\.", "", names(mydat))
prt		<- proc.time()

CL = makeCluster(3)
clusterExport(cl = CL, c("mydat")) 
prt		<- proc.time()
fitls1 	<- parLapply(CL, 1:4, fun = function(cid) {  
require(MCMCglmm)
myprior	<- list(R = list(V = diag(3)/4, n = 3), G = list(G1= list(V = diag(3)/4, n = 3)))
m1		<- MCMCglmm(cbind(y1, y2, y3) ~ trait:X1 + trait:X2 + trait -1, 
					random = ~ idh(trait):id,
					rcov = ~ corg(trait):units, 
					family = c("categorical","categorical","categorical"), 
					prior = myprior, 
					data = mydat,  
					nitt = 11000, burnin = 1000, thin= 40)
})					
use.time2	<- proc.time() - prt
cat("Stan sampling takes", use.time2[3]/60,"min.\n")
stopCluster(CL)

# Diagnosis
beta.mcmc	<- do.call("mcmc.list", lapply(fitls1, function(x) as.mcmc(x$Sol)))
vcv.mcmc	<- do.call("mcmc.list", lapply(fitls1, function(x) as.mcmc(x$VCV)))
plot(beta.mcmc, ask = T)
plot(vcv.mcmc, ask = T)

gelman.diag(beta.mcmc)
gelman.plot(beta.mcmc)
gelman.diag(vcv.mcmc)
gelman.plot(vcv.mcmc)

# Compare with trueth
cbind(HPDinterval(mcmc(as.matrix(beta.mcmc))), postior.mode = posterior.mode(mcmc(as.matrix(beta.mcmc))),true = as.vector(t(beta)))
cbind(HPDinterval(mcmc(as.matrix(vcv.mcmc))), postior.mode = posterior.mode(mcmc(as.matrix(vcv.mcmc))),true = c(omega, Sigma))



