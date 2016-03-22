# Hierarchical MDCEV model 
# 
# Model: 
# U_ht = sum_r psi_{htr}*gamma_r*log(s_{htr}*y_{hr}/(p_{rt}*gamma_r) + 1)
# where 
# log(psi_{htr}) = beta_0{hr} + Xbeta + epsilon_{htr}
# epsilon_{htr} ~ GEV(0, sigma)
# Individual RE: beta_0{hr} ~ N(mu_r, tau_r^2)

library(evd)
library(nloptr)
library(rstan)
library(doParallel)
# source('~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R')
source("0_Allocation_function.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
plot.wd	<- "/kellogg/users/marketing/2661703/Exercise"

#################
# Simulate data # 
#################
# Set simulation parameters
set.seed(66)
nh	<- 20
T	<- 30
N	<- nh*T
R	<- 3
nx	<- 2
beta		<- c(.5, -1)
mu			<- c(0, -.7, .5)
tau			<- c(0, .2, .2)
sigma		<- .7
gamma		<- rep(1, R)
beta0_base	<- 1			# Set the one beta0 = 0
mu[beta0_base]	<- 0
tau[beta0_base]	<- 0

# Simulate data 
idx	<- rep(1:nh, each = T)
X	<- array(runif(R*N*nx, -1, 1), c(R, N, nx))
p	<- matrix(rexp(N*R, 2), N, R)
y	<- rexp(N, .3)
summary(p)
summary(y)
beta0	<- NULL
for(i in 1:R){ beta0	<- cbind(beta0, rnorm(nh, mu[i], tau[i])) }
beta0.ext	<- kronecker(beta0, rep(1, T))
eps		<- matrix(rgev(N*R, scale = sigma), N, R)
psi		<- exp(sapply(1:R, function(i) beta0.ext[,i] + X[i,,]%*%beta + eps[,i]))
e  		<- sapply(1:N, function(i) Allocation_fn(y[i], psi[i,], gamma, Q = Inf, price = p[i,], R, Ra =R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE)$e)
e		<- t(e)
shr		<- e/y
shr_sgn	<- 1*(e > 0)
M		<- rowSums(shr_sgn)

#############
# Stan code # 
#############
src_code	<- '
data{
	int<lower = 1> N; 
	int<lower = 1> nh;
	int<lower = 1> R; 
	int<lower = 1> nx;
	int<lower = 1> beta0_base; 
	int<lower = 1> idx[N]; 
	matrix[N,R]		shr; 
	matrix[N,R]		shr_sgn; 
	matrix[N,R]		p; 
	vector[N]		y; 
	matrix[N,nx]	X[R]; 	
	vector[N]	M; 
	vector[R]	ones; 
	vector[N]	lones;
}
parameters{
	vector[R-1]	mu_raw; 
	vector<lower = 0>[R-1] tau_raw; 
	vector[nx]	beta; 
	vector<lower = 0>[R] gamma; 
	vector[R]	beta0_raw[nh]; 
	real<lower = 0> sigma; 
}
transformed parameters{
	vector[R]	mu; 
	vector[R]	tau; 
	vector[R]	beta0[N]; 
	
	for(j in 1:R){
		if(j==beta0_base){
			mu[j]	<- 0; 
			tau[j]	<- 0;
		}else{
			if(j < beta0_base){
				mu[j]	<- mu_raw[j];
				tau[j]	<- tau_raw[j];
			}else{
				mu[j]	<- mu_raw[(j-1)];
				tau[j]	<- tau_raw[(j-1)]; 
			}
		}
	}
	
	for(i in 1:N){
		for(j in 1:R){
			beta0[i,j]	<- mu[j] + tau[j]*beta0_raw[idx[i],j]; 
		}
	}
}
model{
	// Declare working objects; 
	vector[N]	J; 
	vector[N]	tot; 
	matrix[N,R]		d; 
	matrix[N,R]		V; 
	
	// Prior; 
	mu_raw	~ normal(0, 10);
	tau_raw	~ uniform(0, 10); 
	beta	~ normal(0, 10); 
	gamma	~ uniform(0, 10); 
	for(i in 1:nh){
		beta0_raw[i] ~ normal(0, 1); 
	}
	sigma	~ uniform(0, 10); 
	
	// Likelihood; 
	for(i in 1:N){
		for(j in 1:R){
			d[i,j]	<- shr[i,j] + gamma[j] * p[i,j]/y[i];
			V[i,j]	<- X[j,i]*beta + beta0[i,j] - log(d[i,j]/gamma[j]); 
		}
	}
	J	<- log((d .* shr_sgn)*ones) - (log(d) .* shr_sgn)*ones; 
	tot	<- J - (M-lones)*log(sigma) + (V/sigma .* shr_sgn)*ones -M .* log(exp(V/sigma)*ones); 
	increment_log_prob(sum(tot)); 
}
'

data_ls	<- list(N=N, nh=nh, R=R, nx=nx, beta0_base=beta0_base, idx=idx, shr=shr, shr_sgn=shr_sgn, 
				p=p, y=y, X=X, M=M, ones = rep(1, R), lones = rep(1, N))
save_par<- c("beta", "gamma", "sigma", "mu", "tau")				
init_ls	<- list(mu_raw = as.vector(mu[-beta0_base]), tau_raw = as.vector(tau[-beta0_base]), 
				beta = beta, gamma = gamma, 
				beta_raw = matrix(0, nh, R))
my_init	<- function(){return(init_ls)}
# sim		<- stan(model_code = src_code, data = data_ls, pars = save_par, init = my_init, chains = 4, iter = 5000, thin = 10)

CL = makeCluster(3)
clusterExport(cl = CL, c("src_code", "data_ls", "save_par", "my_init", "init_ls")) 
prt		<- proc.time()
fitls <- parLapply(CL, 1:3, fun = function(cid) {  
  	require(rstan)
	stan(model_code = src_code, data = data_ls, pars = save_par, init = my_init, 
		chains = 1, iter = 100, warmup = 20, thin = 1, chain_id = cid)
})
use.time1	<- proc.time() - prt
stopCluster(CL)
cat("Stan sampling takes", use.time1[3]/60,"min.\n")

fit		<- sflist2stanfit(fitls)
print(fit)
c(beta, gamma, sigma, mu, tau)

