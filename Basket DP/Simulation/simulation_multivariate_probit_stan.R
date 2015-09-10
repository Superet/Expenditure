# Hierarchical multivariate probit

library(mvtnorm)
library(rstan)
library(doParallel)
library(parallel)
library(foreach)
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
y		<- 1*(ystar>0)
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
transformed parameters{
	vector[K]	beta_adj[D]; 
	vector[J]	beta0_adj[D];
	vector[D]	omega_adj;
	corr_matrix[D]	Sigma_adj;
	vector[D]	beta0_mean; 
	
	for(i in 1:D){
		beta0_mean[i]	<- mean(beta0[i]);
		beta_adj[i]		<- beta[i]/sqrt(2);
		beta_adj[i,1]	<- beta_adj[i,1] + beta0_mean[i]/sqrt(2);
		beta0_adj[i]	<- (beta0[i] - beta0_mean[i])/sqrt(2); 
		omega_adj[i]	<- omega[i]/sqrt(2); 
		for(j in 1:D){
			if(i == j){	Sigma_adj[i,j]	<- 1; }
			else{ Sigma_adj[i,j]	<- Sigma[i,j]/sqrt(2); }
		}
	}
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
			y[i,j] ~ bernoulli(Phi(z[i,j])); 
		}
	}
}
"

mydata		<- list(N = N, D = D, K = K, J = J, X = X, y = y, ID = diag(D), idx = id.idx, zz = rep(0, D))
save_par	<- c("beta_adj", "Sigma_adj","omega_adj", "beta0_adj")
# my_init		<- function(){
# 	list(beta = t(as.matrix(beta)), beta0 = t(as.matrix(beta0)), omega = omega, Sigma = Sigma, z = ystar + matrix(rnorm(N*D), N, D) )
# }
init_ls		<- list(beta = t(as.matrix(beta)), beta0 = t(as.matrix(beta0)), omega = omega, Sigma = Sigma, z = ystar + matrix(rnorm(N*D), N, D) )
my_init		<- function(){return(init_ls)}
foo			<- stan(model_code = src_code, data = mydata, pars = save_par, init = my_init, chains = 1, iter = 1)

# Register parallel computing
# mycore 	<- 3
# cl		<- makeCluster(mycore, type = "FORK")
# registerDoParallel(cl)
# cat("Register", mycore, "core parallel computing. \n")
# 
# fitls	<- foreach(i = 1:3) %dopar%{
# 	stan(model_code = src_code, data = mydata, pars = save_par, 
# 						chains = 1, iter = 120, warmup = 20, thin = 4, chain_id = i, refresh = -1 )
# }
# stopCluster(cl)

# fitls <- mclapply(1:3, mc.cores = 3, 
#            function(i) stan(model_code = src_code, data = mydata, pars = save_par, 
# 		 						chains = 1, iter = 120, warmup = 20, thin = 4, chain_id = i, refresh = -1 ))

CL = makeCluster(3)
clusterExport(cl = CL, c("foo", "src_code", "mydata", "save_par", "my_init", "init_ls")) 
fitls <- parLapply(CL, 1:3, fun = function(cid) {  
  	require(rstan)
	stan(fit = foo, data = mydata, pars = save_par, init = my_init,
					chains = 1, iter = 120, warmup = 20, thin = 4, chain_id = cid)
})
stopCluster(CL)

fit		<- sflist2stanfit(fitls)
# 
# fit 		<-  stan(model_code = src_code, data = mydata, pars = save_par, init = my_init,
# 					chains = 3, iter = 120, warmup = 20, thin = 4, cores = getOption("mc.cores",3))

print(fit)

# pdf(paste(plot.wd, "/graph_stanfit.pdf",sep=""), width = 10, height = 8)
# plot(fit)
# dev.off()
